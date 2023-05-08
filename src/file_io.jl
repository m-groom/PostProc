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
    # Initialise viscosity
    μ = zeros(2) # Assumes two species
    report("Reading settings from file $filename")
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
    # Read viscosities
    if (nVars >= 10)
        μ[1] = parse(Float64, split(readline(file), "#")[1])
        μ[2] = parse(Float64, split(readline(file), "#")[1])
    end
    # Close file
    close(file)
    # Package grid settings into a struct
    grid = rectilinearGrid(Nx, Ny, Nz, xL, xR, yL, yR, zL, zR)
    # Package input file settings into a struct
    input = inputSettings(nVars, startTime, Δt, nFiles)

    return grid, input, x0, μ

end

# Function for writing the solution to a .vtr file
function writeSolution(t::Float64, x::Array{Float32,3}, y::Array{Float32,3}, z::Array{Float32,3}, Q::Array{Float32,4}, iNVars::Int64)
    tStart = report("Writing the solution at time $(rpad(string(round(t, digits=5)), 7, "0"))", 1)
    # Check that there are 8, 10 or 12 variables
    if (iNVars == 8) || (iNVars == 10) || (iNVars == 12)
        # Write solution to a .vtr file
        filename = "data/solution_$(rpad(string(round(t, digits=5)), 7, "0")).vtr"
        WriteVTK.vtk_grid(filename, x[:, 1, 1], y[1, :, 1], z[1, 1, :]) do vtk
            vtk["MomentumX", WriteVTK.VTKCellData()] = Q[1:end-1, 1:end-1, 1:end-1, 1]
            vtk["MomentumY", WriteVTK.VTKCellData()] = Q[1:end-1, 1:end-1, 1:end-1, 2]
            vtk["MomentumZ", WriteVTK.VTKCellData()] = Q[1:end-1, 1:end-1, 1:end-1, 3]
            vtk["EnergyDensity", WriteVTK.VTKCellData()] = Q[1:end-1, 1:end-1, 1:end-1, 4]
            if (iNVars >= 10)
                for v = 1:2
                    vtk["DensityMassFraction$(v)", WriteVTK.VTKCellData()] = Q[1:end-1, 1:end-1, 1:end-1, 4+v]
                end
            end
            if (iNVars == 12)
                for v = 1:2
                    vtk["VolumeFraction$(v)", WriteVTK.VTKCellData()] = Q[1:end-1, 1:end-1, 1:end-1, 4+2+v]
                end
            end
            vtk["Density", WriteVTK.VTKCellData()] = Q[1:end-1, 1:end-1, 1:end-1, iNVars-3]
            vtk["Gamma", WriteVTK.VTKCellData()] = Q[1:end-1, 1:end-1, 1:end-1, iNVars-2]
            vtk["Pressure", WriteVTK.VTKCellData()] = Q[1:end-1, 1:end-1, 1:end-1, iNVars-1]
            vtk["Temperature", WriteVTK.VTKCellData()] = Q[1:end-1, 1:end-1, 1:end-1, iNVars]
        end
    else
        error("Number of variables must be 8, 10 or 12")
    end
    tEnd = report("Finished writing solution...", 1)
    report("Elapsed time: $(tEnd - tStart)")
end

# Function for writing out a slice to a .vtr file
function writeSlice(t::Float64, x::Array{Float32,3}, y::Array{Float32,3}, z::Array{Float32,3}, Q::Array{Float32,4}, iNVars::Int64, slice::String, x0::Float64)
    # Calculate y0 and z0
    y0 = (y[1, end, 1] + y[1, 1, 1]) * 0.5
    z0 = (z[1, 1, end] + z[1, 1, 1]) * 0.5
    # Find index of x0, y0 and z0
    iX0 = argmin(abs.(x[:, 1, 1] .- x0))
    iY0 = argmin(abs.(y[1, :, 1] .- y0))
    iZ0 = argmin(abs.(z[1, 1, :] .- z0))
    # Get slice from the solution Q
    if (slice == "xz" || slice == "zx")
        Qs = Q[:, iY0, :, :]
        x1 = x[:, iY0, 1]
        x2 = z[1, iY0, :]
    elseif (slice == "xy" || slice == "yx")
        Qs = Q[:, :, iZ0, :]
        x1 = x[:, 1, iZ0]
        x2 = y[1, :, iZ0]
    elseif (slice == "yz" || slice == "zy")
        Qs = Q[iX0, :, :, :]
        x1 = y[iX0, :, 1]
        x2 = z[iX0, 1, :]
    else
        error("Slice must be 'xy', 'xz' or 'yz'")
    end
    tStart = report("Writing $(slice) slice at time $(rpad(string(round(t, digits=5)), 7, "0"))", 1)
    # Check that there are 8, 10 or 12 variables
    if (iNVars == 8) || (iNVars == 10) || (iNVars == 12)
        # Write slice to a .vtr file
        filename = "data/slice$(uppercase(slice))_$(rpad(string(round(t, digits=5)), 7, "0")).vtr"     
        WriteVTK.vtk_grid(filename, x1, x2) do vtk
            vtk["MomentumX", WriteVTK.VTKCellData()] = Qs[1:end-1, 1:end-1, 1]
            vtk["MomentumY", WriteVTK.VTKCellData()] = Qs[1:end-1, 1:end-1, 2]
            vtk["MomentumZ", WriteVTK.VTKCellData()] = Qs[1:end-1, 1:end-1, 3]
            vtk["EnergyDensity", WriteVTK.VTKCellData()] = Qs[1:end-1, 1:end-1, 4]
            if (iNVars >= 10)
                for v = 1:2
                    vtk["DensityMassFraction$(v)", WriteVTK.VTKCellData()] = Qs[1:end-1, 1:end-1, 4+v]
                end
            end
            if (iNVars == 12)
                for v = 1:2
                    vtk["VolumeFraction$(v)", WriteVTK.VTKCellData()] = Qs[1:end-1, 1:end-1, 4+2+v]
                end
            end
            vtk["Density", WriteVTK.VTKCellData()] = Qs[1:end-1, 1:end-1, iNVars-3]
            vtk["Gamma", WriteVTK.VTKCellData()] = Qs[1:end-1, 1:end-1, iNVars-2]
            vtk["Pressure", WriteVTK.VTKCellData()] = Qs[1:end-1, 1:end-1, iNVars-1]
            vtk["Temperature", WriteVTK.VTKCellData()] = Qs[1:end-1, 1:end-1, iNVars]
        end
    else
        error("Number of variables must be 8, 10 or 12")
    end
    tEnd = report("Finished writing $(slice) slice...", 1)
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
    report("Reading grid partition data from file data/prempi.dat")
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

# Function to write out plane averages to a space delimited text file
function writePlaneAverages(t::Float64, QBar::planeAverage)
    # Make file name
    filename = "data/planeAverages_$(rpad(string(round(t, digits=5)), 7, "0")).dat"
    tStart = report("Writing plane averages to file $filename", 1)
    # Open file
    f = open(filename, "w")
    # Write header
    write(f, "# x   rhoBar   UBar   Y1Bar   Z1Bar   Z1Z2Bar\n")
    # Write data in scientific format with 15 digits
    for i = 1:length(QBar.x)-1
        write(f, "$(@sprintf("%.15e", QBar.x[i]))   $(@sprintf("%.15e", QBar.rhoBar[i]))   $(@sprintf("%.15e", QBar.UBar[i]))   $(@sprintf("%.15e", QBar.Y1Bar[i]))   $(@sprintf("%.15e", QBar.Z1Bar[i]))   $(@sprintf("%.15e", QBar.Z1Z2Bar[i]))\n")
    end
    # Close file
    close(f)
    tEnd = report("Finished writing plane averages", 1)
    report("Elapsed time: $(tEnd - tStart)")
end

# Function to write correlation lengths to file
function writeCorrelationLengths(t::Float64, Λx::Float64, Λyz::Float64)
    # Get filename
    filename = "data/correlationLengths.dat"
    report("Writing correlation lengths to file $filename")
    # Append to file
    f = open(filename, "a")
    # Write header if file is empty
    if (filesize(filename) == 0)
        write(f, "# t   Lambdax   Lambdayz\n")
    end
    # Write data in scientific format with 15 digits
    write(f, "$(@sprintf("%.15e", t))   $(@sprintf("%.15e", Λx))   $(@sprintf("%.15e", Λyz))\n")
    # Close file
    close(f)
end

# Function to write length scales to a space delimited text file
function writeLengthScales(t::Float64, λx::Float64, λyz::Float64, ηx::Float64, ηyz::Float64)
    # Get filename
    filename = "data/lengthScales.dat"
    report("Writing length scales to file $filename")
    # Append to file
    f = open(filename, "a")
    # Write header if file is empty
    if (filesize(filename) == 0)
        write(f, "# t   lambdax   lambdayz   etax   etayz\n")
    end
    # Write data in scientific format with 15 digits
    write(f, "$(@sprintf("%.15e", t))   $(@sprintf("%.15e", λx))   $(@sprintf("%.15e", λyz))   $(@sprintf("%.15e", ηx))   $(@sprintf("%.15e", ηyz))\n")
    # Close file
    close(f)
end

# Function to write out Reynolds stresses to a space delimited text file
function writeReynoldsStresses(t::Float64, x::Array{Float32,1}, R11::Array{Float64,1}, R22::Array{Float64,1}, R33::Array{Float64,1})
    # Make file name
    filename = "data/reynoldsStresses_$(rpad(string(round(t, digits=5)), 7, "0")).dat"
    report("Writing Reynolds stresses to file $filename")
    # Open file
    f = open(filename, "w")
    # Write header
    write(f, "# x   R11   R22   R33\n")
    # Write data in scientific format with 15 digits
    for i = 1:length(x)-1
        write(f, "$(@sprintf("%.15e", x[i]))   $(@sprintf("%.15e", R11[i]))   $(@sprintf("%.15e", R22[i]))   $(@sprintf("%.15e", R33[i]))\n")
    end
    # Close file
    close(f)
end

# Function to write out dissipation rates to a space delimited text file
function writeDissipationRates(t::Float64, x::Array{Float32,1}, εx::Array{Float64,1}, εy::Array{Float64,1}, εz::Array{Float64,1})
    # Make file name
    filename = "data/dissipationRates_$(rpad(string(round(t, digits=5)), 7, "0")).dat"
    report("Writing dissipation rates to file $filename")
    # Open file
    f = open(filename, "w")
    # Write header
    write(f, "# x   epsilonx   epsilony   epsilonz\n")
    # Write data in scientific format with 15 digits
    for i = 1:length(x)-1
        write(f, "$(@sprintf("%.15e", x[i]))   $(@sprintf("%.15e", εx[i]))   $(@sprintf("%.15e", εy[i]))   $(@sprintf("%.15e", εz[i]))\n")
    end
    # Close file
    close(f)
end

# Function to write out Taylor microscales to a space delimited text file
function writeTaylorMicroscales(t::Float64, x::Array{Float32,1}, λx::Array{Float64,1}, λy::Array{Float64,1}, λz::Array{Float64,1})
    # Make file name
    filename = "data/taylor_$(rpad(string(round(t, digits=5)), 7, "0")).dat"
    report("Writing Taylor microscales to file $filename")
    # Open file
    f = open(filename, "w")
    # Write header
    write(f, "# x   lambdax   lambday   lambday\n")
    # Write data in scientific format with 15 digits
    for i = 1:length(x)-1
        write(f, "$(@sprintf("%.15e", x[i]))   $(@sprintf("%.15e", λx[i]))   $(@sprintf("%.15e", λy[i]))   $(@sprintf("%.15e", λz[i]))\n")
    end
    # Close file
    close(f)
end

# Function to write out Kolmogorov microscales to a space delimited text file
function writeKolmogorovMicroscales(t::Float64, x::Array{Float32,1}, ηx::Array{Float64,1}, ηy::Array{Float64,1}, ηz::Array{Float64,1})
    # Make file name
    filename = "data/kolmogorov_$(rpad(string(round(t, digits=5)), 7, "0")).dat"
    report("Writing Kolmogorov microscales to file $filename")
    # Open file
    f = open(filename, "w")
    # Write header
    write(f, "# x   etax   etay   etaz\n")
    # Write data in scientific format with 15 digits
    for i = 1:length(x)-1
        write(f, "$(@sprintf("%.15e", x[i]))   $(@sprintf("%.15e", ηx[i]))   $(@sprintf("%.15e", ηy[i]))   $(@sprintf("%.15e", ηz[i]))\n")
    end
    # Close file
    close(f)
end