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