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
    W = zeros(2) # Assumes two species
    report("Reading settings from file $filename")
    # Open file
    file = open(filename, "r")
    # Read data directory
    dataDir = String(strip(split(readline(file), "#")[1]))
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
    # Read viscosities and molecular weights
    if (nVars >= 10)
        μ[1] = parse(Float64, split(readline(file), "#")[1])
        μ[2] = parse(Float64, split(readline(file), "#")[1])
        W[1] = parse(Float64, split(readline(file), "#")[1])
        W[2] = parse(Float64, split(readline(file), "#")[1])
    end
    # Close file
    close(file)
    # Package grid settings into a struct
    grid = rectilinearGrid(Nx, Ny, Nz, xL, xR, yL, yR, zL, zR)
    # Package input file settings into a struct
    input = inputSettings(nVars, startTime, Δt, nFiles)
    # Package thermodynamic properties into a struct
    thermo = thermodynamicProperties(μ, W)

    return grid, input, thermo, x0, dataDir

end

# Function for writing the solution to a .vtr file
function writeSolution(t::Float64, x::SubArray{Float32,3}, y::SubArray{Float32,3}, z::SubArray{Float32,3}, Q::Array{Float32,4}, iNVars::Int64, dataDir::String)
    tStart = report("Writing the solution at time $(rpad(string(round(t, digits=5)), 7, "0"))", 1)
    # Check that there are 8, 10 or 12 variables
    if (iNVars == 8) || (iNVars == 10) || (iNVars == 12)
        # Write solution to a .vtr file
        filename = "$(dataDir)/solution_$(rpad(string(round(t, digits=5)), 7, "0")).vtr"
        @inbounds WriteVTK.vtk_grid(filename, x[:, 1, 1], y[1, :, 1], z[1, 1, :]) do vtk
            vtk["VelocityX", WriteVTK.VTKCellData()] = @view(Q[1:end-1, 1:end-1, 1:end-1, 1])
            vtk["VelocityY", WriteVTK.VTKCellData()] = @view(Q[1:end-1, 1:end-1, 1:end-1, 2])
            vtk["VelocityZ", WriteVTK.VTKCellData()] = @view(Q[1:end-1, 1:end-1, 1:end-1, 3])
            vtk["EnergyDensity", WriteVTK.VTKCellData()] = @view(Q[1:end-1, 1:end-1, 1:end-1, 4])
            if (iNVars >= 10)
                @inbounds for v = 1:2
                    vtk["MassFraction$(v)", WriteVTK.VTKCellData()] = @view(Q[1:end-1, 1:end-1, 1:end-1, 4+v])
                end
            end
            if (iNVars == 12)
                @inbounds for v = 1:2
                    vtk["VolumeFraction$(v)", WriteVTK.VTKCellData()] = @view(Q[1:end-1, 1:end-1, 1:end-1, 6+v])
                end
            end
            vtk["Density", WriteVTK.VTKCellData()] = @view(Q[1:end-1, 1:end-1, 1:end-1, iNVars-3])
            vtk["Gamma", WriteVTK.VTKCellData()] = @view(Q[1:end-1, 1:end-1, 1:end-1, iNVars-2])
            vtk["Pressure", WriteVTK.VTKCellData()] = @view(Q[1:end-1, 1:end-1, 1:end-1, iNVars-1])
            vtk["Temperature", WriteVTK.VTKCellData()] = @view(Q[1:end-1, 1:end-1, 1:end-1, iNVars])
        end
    else
        error("Number of variables must be 8, 10 or 12")
    end
    tEnd = report("Finished writing solution...", 1)
    report("Elapsed time: $(tEnd - tStart)")
end

# Function for writing out a slice to a .vtr file
function writeSlice(t::Float64, x::SubArray{Float32,3}, y::SubArray{Float32,3}, z::SubArray{Float32,3}, Q::Array{Float32,4}, iNVars::Int64, slice::String, x0::Float64, dataDir::String)
    # Calculate y0 and z0
    y0 = (y[1, end, 1] + y[1, 1, 1]) * 0.5
    z0 = (z[1, 1, end] + z[1, 1, 1]) * 0.5
    # Find index of x0, y0 and z0
    iX0 = searchsortedfirst(x[:, 1, 1], x0) - 1
    iY0 = searchsortedfirst(y[1, :, 1], y0)
    iZ0 = searchsortedfirst(z[1, 1, :], z0)
    # Get slice from the solution Q
    if (slice == "xz" || slice == "zx")
        Qs = @view(Q[:, iY0, :, :])
        x1 = x[:, iY0, 1]
        x2 = z[1, iY0, :]
    elseif (slice == "xy" || slice == "yx")
        Qs = @view(Q[:, :, iZ0, :])
        x1 = x[:, 1, iZ0]
        x2 = y[1, :, iZ0]
    elseif (slice == "yz" || slice == "zy")
        Qs = @view(Q[iX0, :, :, :])
        x1 = y[iX0, :, 1]
        x2 = z[iX0, 1, :]
    else
        error("Slice must be 'xy', 'xz' or 'yz'")
    end
    tStart = report("Writing $(slice) slice at time $(rpad(string(round(t, digits=5)), 7, "0"))", 1)
    # Check that there are 8, 10 or 12 variables
    if (iNVars == 8) || (iNVars == 10) || (iNVars == 12)
        # Write slice to a .vtr file
        filename = "$(dataDir)/slice$(uppercase(slice))_$(rpad(string(round(t, digits=5)), 7, "0")).vtr"
        @inbounds WriteVTK.vtk_grid(filename, x1, x2) do vtk
            vtk["VelocityX", WriteVTK.VTKCellData()] = Qs[1:end-1, 1:end-1, 1]
            vtk["VelocityY", WriteVTK.VTKCellData()] = Qs[1:end-1, 1:end-1, 2]
            vtk["VelocityZ", WriteVTK.VTKCellData()] = Qs[1:end-1, 1:end-1, 3]
            vtk["EnergyDensity", WriteVTK.VTKCellData()] = Qs[1:end-1, 1:end-1, 4]
            if (iNVars >= 10)
                @inbounds for v = 1:2
                    vtk["MassFraction$(v)", WriteVTK.VTKCellData()] = Qs[1:end-1, 1:end-1, 4+v]
                end
            end
            if (iNVars == 12)
                @inbounds for v = 1:2
                    vtk["VolumeFraction$(v)", WriteVTK.VTKCellData()] = Qs[1:end-1, 1:end-1, 6+v]
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
function readPlot3DGrid(timeStep::String, Nx::Int64, Ny::Int64, Nz::Int64, dataDir::String)
    # Read prempi.dat file
    report("Reading file $(dataDir)/prempi.dat")
    NPR, extents = readPrempi("$(dataDir)/prempi.dat")
    # Allocate arrays
    iNX = Int32(extents[1, 2, end] - 2) # Number of points in x-direction
    iNY = Int32(extents[2, 2, end] - 2) # Number of points in y-direction
    iNZ = Int32(extents[3, 2, end] - 2) # Number of points in z-direction
    # Check that the grid sizes match
    if (iNX != Nx + 1) || (iNY != Ny + 1) || (iNZ != Nz + 1)
        error("Grid sizes from $(dataDir)/prempi.dat do not match those in post.par")
    end
    data = zeros(Float32, iNX, iNY, iNZ, 3)
    # Read grid file
    tStart = report("Reading grid files for $(NPR) processors located in $(dataDir)/out_$(timeStep)", 1)
    @inbounds for n = 1:NPR
        # Generate filename
        proc = string(n - 1)
        filename = "$(dataDir)/out_$(timeStep)/$(timeStep).$(proc).g"
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
            error("Grid sizes from $(filename) do not match those in $(dataDir)/prempi.dat")
        end
        data[iStart:iEnd, jStart:jEnd, kStart:kEnd, :] = read(f, (Float32, (ie, je, ke, 3)))
        # Close file
        close(f)
    end
    tEnd = report("Finished reading grid files...", 1)
    report("Elapsed time: $(tEnd - tStart)")
    # Extract x, y and z
    x = @view(data[:, :, :, 1])
    y = @view(data[:, :, :, 2])
    z = @view(data[:, :, :, 3])

    return x, y, z

end

# Function for reading a Plot3D solution file
function readPlot3DSolution(timeStep::String, Nx::Int64, Ny::Int64, Nz::Int64, nVars::Int64, dataDir::String)
    # Read prempi.dat file
    report("Reading grid partition data from file $(dataDir)/prempi.dat")
    NPR, extents = readPrempi("$(dataDir)/prempi.dat")
    # Allocate arrays
    iNX = Int32(extents[1, 2, end] - 2) # Number of points in x-direction
    iNY = Int32(extents[2, 2, end] - 2) # Number of points in y-direction
    iNZ = Int32(extents[3, 2, end] - 2) # Number of points in z-direction
    # Check that the grid sizes match
    if (iNX != Nx + 1) || (iNY != Ny + 1) || (iNZ != Nz + 1)
        error("Grid sizes from $(dataDir)/prempi.dat do not match those in post.par")
    end
    Q = zeros(Float32, iNX, iNY, iNZ, nVars)
    # Read solution file
    tStart = report("Reading solution files for $(NPR) processors located in $(dataDir)/out_$(timeStep)", 1)
    @inbounds for n = 1:NPR
        # Generate filename
        proc = string(n - 1)
        filename = "$(dataDir)/out_$(timeStep)/$(timeStep).$(proc).all.f"
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
            error("Grid sizes from $(filename) do not match those in $(dataDir)/prempi.dat")
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
    @inbounds for n = 1:NPR
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
function writePlaneAverages(t::Float64, QBar::planeAverage, grid::rectilinearGrid, dataDir::String)
    # Find index where x = xL
    iL = searchsortedfirst(QBar.x, grid.xL)
    # Find index where x = xR
    iR = searchsortedfirst(QBar.x, grid.xR)
    # Make file name
    filename = "$(dataDir)/planeAverages_$(rpad(string(round(t, digits=5)), 7, "0")).dat"
    tStart = report("Writing plane averages to file $filename", 1)
    # Open file
    f = open(filename, "w")
    # Write header
    write(f, "# x   rhoBar   UBar   Y1Bar   Z1Bar   Z1Z2Bar\n")
    # Write data in scientific format with 15 digits
    @inbounds for i = iL:iR
        write(f, "$(@sprintf("%.15e", QBar.x[i]))   $(@sprintf("%.15e", QBar.rhoBar[i]))   $(@sprintf("%.15e", QBar.UBar[i]))   $(@sprintf("%.15e", QBar.Y1Bar[i]))   $(@sprintf("%.15e", QBar.Z1Bar[i]))   $(@sprintf("%.15e", QBar.Z1Z2Bar[i]))\n")
    end
    # Close file
    close(f)
    tEnd = report("Finished writing plane averages", 1)
    report("Elapsed time: $(tEnd - tStart)")
end

# Function to write correlation lengths to file
function writeCorrelationLengths(t::Float64, x::SubArray{Float32,1}, Λx::Array{Float64,1}, Λyz::Array{Float64,1}, grid::rectilinearGrid, dataDir::String)
    # Find index where x = xL
    iL = searchsortedfirst(x, grid.xL)
    # Find index where x = xR
    iR = searchsortedfirst(x, grid.xR)
    # Make file name
    filename = "$(dataDir)/correlationLengths_$(rpad(string(round(t, digits=5)), 7, "0")).dat"
    report("Writing correlation lengths to file $filename")
    # Open file
    f = open(filename, "w")
    # Write header
    write(f, "# x   Lambdax   Lambdayz\n")
    # Write data in scientific format with 15 digits
    @inbounds for i = iL:iR
        write(f, @sprintf("%.15e", x[i]),"   ", @sprintf("%.15e", Λx[i]),"   ", @sprintf("%.15e", Λyz[i]), "\n")
    end
    # Close file
    close(f)
end

# Function to write length scales to a space delimited text file
function writeLengthScales(t::Float64, λx::Float64, λyz::Float64, ηx::Float64, ηyz::Float64, dataDir::String)
    # Get filename
    filename = "$(dataDir)/lengthScales.dat"
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
function writeReynoldsStresses(t::Float64, x::SubArray{Float32,1}, R11::Array{Float64,1}, R22::Array{Float64,1}, R33::Array{Float64,1}, grid::rectilinearGrid, dataDir::String)
    # Find index where x = xL
    iL = searchsortedfirst(x, grid.xL)
    # Find index where x = xR
    iR = searchsortedfirst(x, grid.xR)
    # Make file name
    filename = "$(dataDir)/reynoldsStresses_$(rpad(string(round(t, digits=5)), 7, "0")).dat"
    report("Writing Reynolds stresses to file $filename")
    # Open file
    f = open(filename, "w")
    # Write header
    write(f, "# x   R11   R22   R33\n")
    # Write data in scientific format with 15 digits
    @inbounds for i = iL:iR
        write(f, "$(@sprintf("%.15e", x[i]))   $(@sprintf("%.15e", R11[i]))   $(@sprintf("%.15e", R22[i]))   $(@sprintf("%.15e", R33[i]))\n")
    end
    # Close file
    close(f)
end

# Function to write out dissipation rates to a space delimited text file
function writeDissipationRates(t::Float64, x::SubArray{Float32,1}, εx::Array{Float64,1}, εy::Array{Float64,1}, εz::Array{Float64,1}, grid::rectilinearGrid, dataDir::String)
    # Find index where x = xL
    iL = searchsortedfirst(x, grid.xL)
    # Find index where x = xR
    iR = searchsortedfirst(x, grid.xR)
    # Make file name
    filename = "$(dataDir)/dissipationRates_$(rpad(string(round(t, digits=5)), 7, "0")).dat"
    report("Writing dissipation rates to file $filename")
    # Open file
    f = open(filename, "w")
    # Write header
    write(f, "# x   epsilonx   epsilony   epsilonz\n")
    # Write data in scientific format with 15 digits
    @inbounds for i = iL:iR
        write(f, "$(@sprintf("%.15e", x[i]))   $(@sprintf("%.15e", εx[i]))   $(@sprintf("%.15e", εy[i]))   $(@sprintf("%.15e", εz[i]))\n")
    end
    # Close file
    close(f)
end

# Function to write out Taylor microscales to a space delimited text file
function writeTaylorMicroscales(t::Float64, x::SubArray{Float32,1}, λx::Array{Float64,1}, λy::Array{Float64,1}, λz::Array{Float64,1}, grid::rectilinearGrid, dataDir::String)
    # Find index where x = xL
    iL = searchsortedfirst(x, grid.xL)
    # Find index where x = xR
    iR = searchsortedfirst(x, grid.xR)
    # Make file name
    filename = "$(dataDir)/taylor_$(rpad(string(round(t, digits=5)), 7, "0")).dat"
    report("Writing Taylor microscales to file $filename")
    # Open file
    f = open(filename, "w")
    # Write header
    write(f, "# x   lambdax   lambday   lambday\n")
    # Write data in scientific format with 15 digits
    @inbounds for i = iL:iR
        write(f, "$(@sprintf("%.15e", x[i]))   $(@sprintf("%.15e", λx[i]))   $(@sprintf("%.15e", λy[i]))   $(@sprintf("%.15e", λz[i]))\n")
    end
    # Close file
    close(f)
end

# Function to write out Kolmogorov microscales to a space delimited text file
function writeKolmogorovMicroscales(t::Float64, x::SubArray{Float32,1}, ηx::Array{Float64,1}, ηy::Array{Float64,1}, ηz::Array{Float64,1}, grid::rectilinearGrid, dataDir::String)
    # Find index where x = xL
    iL = searchsortedfirst(x, grid.xL)
    # Find index where x = xR
    iR = searchsortedfirst(x, grid.xR)
    # Make file name
    filename = "$(dataDir)/kolmogorov_$(rpad(string(round(t, digits=5)), 7, "0")).dat"
    report("Writing Kolmogorov microscales to file $filename")
    # Open file
    f = open(filename, "w")
    # Write header
    write(f, "# x   etax   etay   etaz\n")
    # Write data in scientific format with 15 digits
    @inbounds for i = iL:iR
        write(f, "$(@sprintf("%.15e", x[i]))   $(@sprintf("%.15e", ηx[i]))   $(@sprintf("%.15e", ηy[i]))   $(@sprintf("%.15e", ηz[i]))\n")
    end
    # Close file
    close(f)
end

# Function to write integral width to a space delimited text file
function writeIntegralWidth(t::Float64, W::Float64, H::Float64, dataDir::String)
    # Get filename
    filename = "$(dataDir)/integralWidth.dat"
    report("Writing integral and product widths to file $filename")
    # Append to file
    f = open(filename, "a")
    # Write header if file is empty
    if (filesize(filename) == 0)
        write(f, "# t   W   H\n")
    end
    # Write data in scientific format with 15 digits
    write(f, "$(@sprintf("%.15e", t))   $(@sprintf("%.15e", W))   $(@sprintf("%.15e", H))\n")
    # Close file
    close(f)
end

# Function to write energy spectra to a space delimited text file
function writeEnergySpectra(t::Float64, κ::Array{Float64,1}, Ex::Array{Float64,1}, Ey::Array{Float64,1}, Ez::Array{Float64,1}, grid::rectilinearGrid, num::Array{Int64,1}, dataDir::String)
    # Calculate denoising term
    Δκy = 2.0 * π / (grid.yR - grid.yL)
    Δκz = 2.0 * π / (grid.zR - grid.zL)
    Δκ = min(Δκy, Δκz)
    denoise = 2.0 .* π .* κ ./ (Δκ * num)
    # Denoise spectra
    Ex = Ex .* denoise
    Ey = Ey .* denoise
    Ez = Ez .* denoise
    # Calculate Nyquist wavenumber
    N = Int(max(grid.Ny / 2, grid.Nz / 2))
    # Make file name
    filename = "$(dataDir)/spectra_$(rpad(string(round(t, digits=5)), 7, "0")).dat"
    report("Writing energy spectra to file $filename")
    # Open file
    f = open(filename, "w")
    # Write header
    write(f, "# kappa   Ex   Ey   Ez\n")
    # Write data in scientific format with 15 digits
    @inbounds for i = 1:N+1
        write(f, "$(@sprintf("%.15e", κ[i]))   $(@sprintf("%.15e", Ex[i]))   $(@sprintf("%.15e", Ey[i]))   $(@sprintf("%.15e", Ez[i]))\n")
    end
    # Close file
    close(f)
end

# Function to write integral length to a space delimited text file
function writeIntegralLength(t::Float64, Lyz::Float64, dataDir::String)
    # Get filename
    filename = "$(dataDir)/integralLength.dat"
    report("Writing integral length to file $filename")
    # Append to file
    f = open(filename, "a")
    # Write header if file is empty
    if (filesize(filename) == 0)
        write(f, "# t   Lyz\n")
    end
    # Write data in scientific format with 15 digits
    write(f, "$(@sprintf("%.15e", t))   $(@sprintf("%.15e", Lyz))\n")
    # Close file
    close(f)
end

# Function to write out velocity correlation to a space delimited text file
function writeVelocityCorrelation(t::Float64, x::SubArray{Float32,1}, R11::Array{Float64,2}, R22::Array{Float64,2}, R33::Array{Float64,2}, grid::rectilinearGrid, dataDir::String)
    # Find index where x = xL
    iL = searchsortedfirst(x, grid.xL)
    # Find index where x = xR
    iR = searchsortedfirst(x, grid.xR)
    # Make file name for R11
    filename = "$(dataDir)/velocityCorrelationX_$(rpad(string(round(t, digits=5)), 7, "0")).dat"
    report("Writing velocity correlation R11 to file $filename")
    # Open file
    f1 = open(filename, "w")
    # Write header
    write(f1, "# x   R11\n")
    # Write data in scientific format with 15 digits
    @inbounds for i = iL:iR
        write(f1, @sprintf("%.15e", x[i]), "   ")
        @inbounds for j = 1:size(R11, 1)
            write(f1, @sprintf("%.15e", R11[j, i]), "   ")
        end
        write(f1, "\n")
    end
    # Close file
    close(f1)
    # Make file name for R22
    filename = "$(dataDir)/velocityCorrelationY_$(rpad(string(round(t, digits=5)), 7, "0")).dat"
    report("Writing velocity correlation R22 to file $filename")
    # Open file
    f2 = open(filename, "w")
    # Write header
    write(f2, "# x   R22\n")
    # Write data in scientific format with 15 digits
    @inbounds for i = iL:iR
        write(f2, @sprintf("%.15e", x[i]), "   ")
        @inbounds for j = 1:size(R22, 1)
            write(f2, @sprintf("%.15e", R22[j, i]), "   ")
        end
        write(f2, "\n")
    end
    # Close file
    close(f2)
    # Make file name for R33
    filename = "$(dataDir)/velocityCorrelationZ_$(rpad(string(round(t, digits=5)), 7, "0")).dat"
    report("Writing velocity correlation R33 to file $filename")
    # Open file
    f3 = open(filename, "w")
    # Write header
    write(f3, "# x   R33\n")
    # Write data in scientific format with 15 digits
    @inbounds for i = iL:iR
        write(f3, @sprintf("%.15e", x[i]), "   ")
        @inbounds for j = 1:size(R33, 1)
            write(f3, @sprintf("%.15e", R33[j, i]), "   ")
        end
        write(f3, "\n")
    end
    # Close file
    close(f3)
end