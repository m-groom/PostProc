# Functions to for file input and output

# Function for printing a message to stdout and a log file
function report(message, returnTime=0)
    # Get the current time
    now = Dates.now()
    # Prepend the current time to the message
    message = "$(Dates.format(now, "yyyy-mm-dd HH:MM:SS.sss"))   $(message)"
    # Print the message to stdout
    println(message)
    # Write the message to the log file
    file = open("processing.log", "a")
    println(file, message)
    close(file)
    if (returnTime == 1)
        return now
    end
end

# function to read settings for the post.par file
function readSettings(filename::String)
    # Open file
    file = open(filename, "r")
    # Read grid settings
    Nx = parse(Int64, split(readline(file), "#")[1])
    Ny = parse(Int64, split(readline(file), "#")[1])
    Nz = parse(Int64, split(readline(file), "#")[1])
    xL = parse(Float64, split(readline(file), "#")[1])
    xR = parse(Float64, split(readline(file), "#")[1])
    yL = parse(Float64, split(readline(file), "#")[1])
    yR = parse(Float64, split(readline(file), "#")[1])
    zL = parse(Float64, split(readline(file), "#")[1])
    zR = parse(Float64, split(readline(file), "#")[1])
    x0 = parse(Float64, split(readline(file), "#")[1])
    # Read solution settings
    nVars = parse(Int, split(readline(file), "#")[1])
    startTime = parse(Float64, split(readline(file), "#")[1])
    Δt = parse(Float64, split(readline(file), "#")[1])
    nFiles = parse(Int, split(readline(file), "#")[1])
    # Close file
    close(file)
    # Package grid settings into a struct
    grid = rectilinearGrid(Nx, Ny, Nz, xL, xR, yL, yR, zL, zR)
    # Package input file settings into a struct
    input = inputSettings(nVars, startTime, Δt, nFiles)

    return grid, input, x0

end

# Function for writing the solution to a .vtr file
function writeSolution(t::Float64, x::Array{Float32,3}, y::Array{Float32,3}, z::Array{Float32,3}, Q::Array{Float32,4}, iNVars::Int64)
    # Check that there are 8, 10 or 12 variables
    if (iNVars == 8) || (iNVars == 10) || (iNVars == 12)
        # Write solution to a vtr file
        filename = "data/solution_$(rpad(string(round(t, digits=5)), 7, "0")).vtr"
        tStart = report("Writing the solution at time $(rpad(string(round(t, digits=5)), 7, "0"))", 1)
        WriteVTK.vtk_grid(filename, x[:, 1, 1], y[1, :, 1], z[1, 1, :]) do vtk
            vtk["MomentumX"] = Q[:, :, :, 1]
            vtk["MomentumY"] = Q[:, :, :, 2]
            vtk["MomentumZ"] = Q[:, :, :, 3]
            vtk["EnergyDensity"] = Q[:, :, :, 4]
            if (iNVars >= 10)
                for v = 1:2
                    vtk["DensityMassFraction$(v)"] = Q[:, :, :, 4+v]
                end
            end
            if (iNVars == 12)
                for v = 1:2
                    vtk["VolumeFraction$(v)"] = Q[:, :, :, 4+2+v]
                end
            end
            vtk["Density"] = Q[:, :, :, iNVars-3]
            vtk["Gamma"] = Q[:, :, :, iNVars-2]
            vtk["Pressure"] = Q[:, :, :, iNVars-1]
            vtk["Temperature"] = Q[:, :, :, iNVars]
        end
    else
        error("Number of variables must be 8, 10 or 12")
    end
    tEnd = report("Finished writing solution...", 1)
    report("Elapsed time: $(tEnd - tStart)")
end

# Function for reading a Plot3D grid file
function readPlot3DGrid(timeStep::String, Nx::Int64, Ny::Int64, Nz::Int64)
    # Read prempi.dat file
    report("Reading file data/prempi.dat")
    NPR, extents = readPrempi("data/prempi.dat")
    # Allocate arrays
    iNX = Int32(extents[1, 2, end] - 2) # Number of points in x-direction
    iNY = Int32(extents[2, 2, end] - 2) # Number of points in y-direction
    iNZ = Int32(extents[3, 2, end] - 2) # Number of points in z-direction
    # Check that the grid sizes match
    if (iNX != Nx + 1) || (iNY != Ny + 1) || (iNZ != Nz + 1)
        error("Grid sizes from data/prempi.dat do not match those in post.par")
    end
    data = Float32.(zeros(iNX, iNY, iNZ, 3))
    # Read grid file
    tStart = report("Reading grid files for $(NPR) processors located in data/out_$(timeStep)", 1)
    for n = 1:NPR
        # Generate filename
        proc = string(n - 1)
        filename = "data/out_$(timeStep)/$(timeStep).$(proc).g"
        # Get extents for this processor
        iStart = Int32(extents[1, 1, n])
        iEnd = Int32(extents[1, 2, n] - 2)
        jStart = Int32(extents[2, 1, n])
        jEnd = Int32(extents[2, 2, n] - 2)
        kStart = Int32(extents[3, 1, n])
        kEnd = Int32(extents[3, 2, n] - 2)
        # Open file
        report("Reading grid file for processor $(proc)")
        f = FFile.FortranFile(filename, "r")
        blocks = read(f, Int32)
        ie, je, ke = read(f, (Int32, 3))
        # Check that ie, je, ke match
        if (ie != iEnd - iStart + 1) || (je != jEnd - jStart + 1) || (ke != kEnd - kStart + 1)
            error("Grid sizes from $(filename) do not match those in data/prempi.dat")
        end
        data[iStart:iEnd, jStart:jEnd, kStart:kEnd, :] = read(f, (Float32, (ie, je, ke, 3)))
        # Close file
        close(f)
    end
    tEnd = report("Finished reading grid files...", 1)
    report("Elapsed time: $(tEnd - tStart)")

    return data[:, :, :, 1], data[:, :, :, 2], data[:, :, :, 3]

end

# Function for reading a Plot3D solution file
function readPlot3DSolution(timeStep::String, Nx::Int64, Ny::Int64, Nz::Int64, nVars::Int64)
    # Read prempi.dat file
    report("Reading file data/prempi.dat")
    NPR, extents = readPrempi("data/prempi.dat")
    # Allocate arrays
    iNX = Int32(extents[1, 2, end] - 2) # Number of points in x-direction
    iNY = Int32(extents[2, 2, end] - 2) # Number of points in y-direction
    iNZ = Int32(extents[3, 2, end] - 2) # Number of points in z-direction
    # Check that the grid sizes match
    if (iNX != Nx + 1) || (iNY != Ny + 1) || (iNZ != Nz + 1)
        error("Grid sizes from data/prempi.dat do not match those in post.par")
    end
    Q = Float32.(zeros(iNX, iNY, iNZ, nVars))
    # Read solution file
    tStart = report("Reading solution files for $(NPR) processors located in data/out_$(timeStep)", 1)
    for n = 1:NPR
        # Generate filename
        proc = string(n - 1)
        filename = "data/out_$(timeStep)/$(timeStep).$(proc).all.f"
        # Get extents for this processor
        iStart = Int32(extents[1, 1, n])
        iEnd = Int32(extents[1, 2, n] - 2)
        jStart = Int32(extents[2, 1, n])
        jEnd = Int32(extents[2, 2, n] - 2)
        kStart = Int32(extents[3, 1, n])
        kEnd = Int32(extents[3, 2, n] - 2)
        # Open file
        report("Reading solution file for processor $(proc)")
        f = FFile.FortranFile(filename, "r")
        blocks = read(f, Int32)
        ie, je, ke = read(f, (Int32, 3))
        # Check that ie, je, ke match
        if (ie != iEnd - iStart + 1) || (je != jEnd - jStart + 1) || (ke != kEnd - kStart + 1)
            error("Grid sizes from $(filename) do not match those in data/prempi.dat")
        end
        Q[iStart:iEnd, jStart:jEnd, kStart:kEnd, :] = read(f, (Float32, (ie, je, ke, nVars)))
        # Close file
        close(f)
    end
    tEnd = report("Finished reading solution files...", 1)
    report("Elapsed time: $(tEnd - tStart)")

    return Q

end

# Function for reading a FLAMENCO prempi.dat file
function readPrempi(filename::String)
    # Open file
    file = open(filename, "r")
    # Read number of processors
    NPR = parse(Int, readline(file))
    extents = Int.(zeros(3, 2, NPR))
    for n = 1:NPR
        # Read processor number
        proc = parse(Int, readline(file))
        # Read k extents
        extents[3, :, proc+1] = parse.(Int, split(readline(file)))
        # Read j extents
        extents[2, :, proc+1] = parse.(Int, split(readline(file)))
        # Read i extents
        extents[1, :, proc+1] = parse.(Int, split(readline(file)))
    end
    # Close file
    close(file)

    return NPR, extents

end