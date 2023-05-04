# Functions to read and write VTK files

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
    # Read solution settings
    nVars = parse(Int, split(readline(file), "#")[1])
    startTime = parse(Float64, split(readline(file), "#")[1])
    ΔT = parse(Float64, split(readline(file), "#")[1])
    nFiles = parse(Int, split(readline(file), "#")[1])
    # Close file
    close(file)
    # Package grid settings into a struct
    grid = rectilinearGrid(Nx, Ny, Nz, xL, xR, yL, yR, zL, zR)
    # Package input file settings into a struct
    input = inputSettings(nVars, startTime, ΔT, nFiles)

    return grid, input

end

# Function for reading a Plot3D grid file
function readPlot3DGrid(timeStep::String, Nx::Int64, Ny::Int64, Nz::Int64)
    # Read prempi.dat file
    NPR, extents = readPrempi("data/prempi.dat")
    # Allocate arrays
    iNX = Int32(extents[1, 2, end] - 2) # Number of points in x-direction
    iNY = Int32(extents[2, 2, end] - 2) # Number of points in y-direction
    iNZ = Int32(extents[3, 2, end] - 2) # Number of points in z-direction
    # Check that the grid sizes match
    if (iNX != Nx + 1) || (iNY != Ny + 1) || (iNZ != Nz + 1)
        error("Grid sizes from data/prempi.dat do not match those in grid.par")
    end
    x = Float32.(zeros(iNX, iNY, iNZ))
    y = Float32.(zeros(iNX, iNY, iNZ))
    z = Float32.(zeros(iNX, iNY, iNZ))
    # Read grid file
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
        println("Reading grid file for processor ", proc)
        f = FFile.FortranFile(filename, "r")
        blocks = read(f, Int32)
        ie, je, ke = read(f, (Int32, 3))
        # Check that ie, je, ke match
        if (ie != iEnd-iStart+1) || (je != jEnd-jStart+1) || (ke != kEnd-kStart+1)
            error("Grid sizes from $(filename) do not match those in data/prempi.dat")
        end
        data = read(f, (Float32, (ie, je, ke, 3)))
        x[iStart:iEnd, jStart:jEnd, kStart:kEnd] = data[:, :, :, 1]
        y[iStart:iEnd, jStart:jEnd, kStart:kEnd] = data[:, :, :, 2]
        z[iStart:iEnd, jStart:jEnd, kStart:kEnd] = data[:, :, :, 3]
        # Close file
        close(f)
    end

    return x, y, z

end

# Function for reading a Plot3D solution file
function readPlot3DSolution(timeStep::String, Nx::Int64, Ny::Int64, Nz::Int64, nVars::Int64)
    # Read prempi.dat file
    NPR, extents = readPrempi("data/prempi.dat")
    # Allocate arrays
    iNX = Int32(extents[1, 2, end] - 2) # Number of points in x-direction
    iNY = Int32(extents[2, 2, end] - 2) # Number of points in y-direction
    iNZ = Int32(extents[3, 2, end] - 2) # Number of points in z-direction
    # Check that the grid sizes match
    if (iNX != Nx + 1) || (iNY != Ny + 1) || (iNZ != Nz + 1)
        error("Grid sizes from data/prempi.dat do not match those in grid.par")
    end
    Q = Float32.(zeros(iNX, iNY, iNZ, nVars))
    # Read grid file
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
        println("Reading solution file for processor ", proc, " at time ", timeStep)
        f = FFile.FortranFile(filename, "r")
        blocks = read(f, Int32)
        ie, je, ke = read(f, (Int32, 3))
        # Check that ie, je, ke match
        if (ie != iEnd-iStart+1) || (je != jEnd-jStart+1) || (ke != kEnd-kStart+1)
            error("Grid sizes from $(filename) do not match those in data/prempi.dat")
        end
        Q[iStart:iEnd, jStart:jEnd, kStart:kEnd, :] = read(f, (Float32, (ie, je, ke, nVars)))
        # Close file
        close(f)
    end

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

# Struct to store settings for a rectilinear grid
struct rectilinearGrid
    Nx::Int64
    Ny::Int64
    Nz::Int64
    xL::Float64
    xR::Float64
    yL::Float64
    yR::Float64
    zL::Float64
    zR::Float64
end

# Struct to store input file settings
struct inputSettings
    nVars::Int
    startTime::Float64
    ΔT::Float64
    nFiles::Int
end