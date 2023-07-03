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
function readSettings(filename::AbstractString)
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
function writeSolution(t::AbstractFloat, x::AbstractArray, y::AbstractArray, z::AbstractArray, Q::AbstractArray, iNVars::Integer, dataDir::AbstractString)
    tStart = report("Writing the solution at time $(rpad(string(round(t, digits=5)), 7, "0"))", 1)
    # Check that there are 8, 10 or 12 variables
    if (iNVars == 8) || (iNVars == 10) || (iNVars == 12)
        # Write solution to a .vtr file
        filename = "$(dataDir)/solution_$(rpad(string(round(t, digits=5)), 7, "0")).vtr"
        @inbounds begin
            WriteVTK.vtk_grid(filename, x[:, 1, 1], y[1, :, 1], z[1, 1, :]) do vtk
                vtk["VelocityX", WriteVTK.VTKCellData()] = @view(Q[1:end-1, 1:end-1, 1:end-1, 1])
                vtk["VelocityY", WriteVTK.VTKCellData()] = @view(Q[1:end-1, 1:end-1, 1:end-1, 2])
                vtk["VelocityZ", WriteVTK.VTKCellData()] = @view(Q[1:end-1, 1:end-1, 1:end-1, 3])
                if (iNVars >= 10)
                    vtk["MassFraction1", WriteVTK.VTKCellData()] = @view(Q[1:end-1, 1:end-1, 1:end-1, 5])
                end
                if (iNVars == 12)
                    vtk["VolumeFraction1", WriteVTK.VTKCellData()] = @view(Q[1:end-1, 1:end-1, 1:end-1, 7])
                end
                vtk["Density", WriteVTK.VTKCellData()] = @view(Q[1:end-1, 1:end-1, 1:end-1, iNVars-3])
                vtk["Pressure", WriteVTK.VTKCellData()] = @view(Q[1:end-1, 1:end-1, 1:end-1, iNVars-1])
                vtk["Temperature", WriteVTK.VTKCellData()] = @view(Q[1:end-1, 1:end-1, 1:end-1, iNVars])
            end
        end
    else
        error("Number of variables must be 8, 10 or 12")
    end
    tEnd = report("Finished writing solution...", 1)
    report("Elapsed time: $(tEnd - tStart)")
end

# Function for writing out a slice to a .vtr file
function writeSlice(t::AbstractFloat, x::AbstractArray, y::AbstractArray, z::AbstractArray, Q::AbstractArray, iNVars::Int64, slice::AbstractString, x0::AbstractFloat, dataDir::AbstractString)
    # Calculate y0 and z0
    y0 = (y[1, end, 1] + y[1, 1, 1]) * 0.5
    z0 = (z[1, 1, end] + z[1, 1, 1]) * 0.5
    # Find index of x0, y0 and z0
    iX0 = searchsortedfirst(x[:, 1, 1], x0)
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
        @inbounds begin
            WriteVTK.vtk_grid(filename, x1, x2) do vtk
                vtk["VelocityX", WriteVTK.VTKCellData()] = Qs[1:end-1, 1:end-1, 1]
                vtk["VelocityY", WriteVTK.VTKCellData()] = Qs[1:end-1, 1:end-1, 2]
                vtk["VelocityZ", WriteVTK.VTKCellData()] = Qs[1:end-1, 1:end-1, 3]
                if (iNVars >= 10)
                    vtk["MassFraction1", WriteVTK.VTKCellData()] = Qs[1:end-1, 1:end-1, 5]
                end
                if (iNVars == 12)
                    vtk["VolumeFraction1", WriteVTK.VTKCellData()] = Qs[1:end-1, 1:end-1, 7]
                end
                vtk["Density", WriteVTK.VTKCellData()] = Qs[1:end-1, 1:end-1, iNVars-3]
                vtk["Pressure", WriteVTK.VTKCellData()] = Qs[1:end-1, 1:end-1, iNVars-1]
                vtk["Temperature", WriteVTK.VTKCellData()] = Qs[1:end-1, 1:end-1, iNVars]
            end
        end
    else
        error("Number of variables must be 8, 10 or 12")
    end
    tEnd = report("Finished writing $(slice) slice...", 1)
    report("Elapsed time: $(tEnd - tStart)")
end

# Function for reading a Plot3D grid file
function readPlot3DGrid(timeStep::AbstractString, Nx::Int64, Ny::Int64, Nz::Int64, dataDir::AbstractString)
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
    @inbounds begin
        @threads for n = 1:NPR
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
function readPlot3DSolution(timeStep::AbstractString, Nx::Int64, Ny::Int64, Nz::Int64, nVars::Int64, dataDir::AbstractString)
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
    @inbounds begin
        @threads for n = 1:NPR
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
    end
    tEnd = report("Finished reading solution files...", 1)
    report("Elapsed time: $(tEnd - tStart)")

    return Q

end

# Function for reading a FLAMENCO prempi.dat file
function readPrempi(filename::AbstractString)
    # Open file
    file = open(filename, "r")
    # Read number of processors
    NPR = parse(Int, readline(file))
    extents = Int.(zeros(3, 2, NPR))
    @inbounds begin
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
    end
    # Close file
    close(file)

    return NPR, extents

end

# Function to write out plane averages to a space delimited text file
function writePlaneAverages(t::AbstractFloat, QBar::planeAverage, grid::rectilinearGrid, dataDir::AbstractString)
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
    @inbounds begin
        for i = iL:iR
            write(f, @sprintf("%.15e", QBar.x[i]), "   ", @sprintf("%.15e", QBar.rhoBar[i]), "   ", @sprintf("%.15e", QBar.UBar[i]), "   ", @sprintf("%.15e", QBar.Y1Bar[i]), "   ", @sprintf("%.15e", QBar.Z1Bar[i]), "   ", @sprintf("%.15e", QBar.Z1Z2Bar[i]), "\n")
        end
    end
    # Close file
    close(f)
    tEnd = report("Finished writing plane averages", 1)
    report("Elapsed time: $(tEnd - tStart)")
end

# Function to write out velocity correlation to a space delimited text file
function writeVelocityCorrelation(t::AbstractFloat, x::AbstractArray, R11::AbstractArray, R22::AbstractArray, R33::AbstractArray, grid::rectilinearGrid, dataDir::AbstractString)
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
    @inbounds begin
        for i = iL:iR
            write(f1, @sprintf("%.15e", x[i]), "   ")
            for j = 1:size(R11, 1)
                write(f1, @sprintf("%.15e", R11[j, i]), "   ")
            end
            write(f1, "\n")
        end
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
    @inbounds begin
        for i = iL:iR
            write(f2, @sprintf("%.15e", x[i]), "   ")
            for j = 1:size(R22, 1)
                write(f2, @sprintf("%.15e", R22[j, i]), "   ")
            end
            write(f2, "\n")
        end
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
    @inbounds begin
        for i = iL:iR
            write(f3, @sprintf("%.15e", x[i]), "   ")
            for j = 1:size(R33, 1)
                write(f3, @sprintf("%.15e", R33[j, i]), "   ")
            end
            write(f3, "\n")
        end
    end
    # Close file
    close(f3)
end

# Function to write correlation lengths to file
function writeCorrelationLengths(t::AbstractFloat, x::AbstractArray, Λx::AbstractArray, Λy::AbstractArray, Λz::AbstractArray, grid::rectilinearGrid, dataDir::AbstractString)
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
    write(f, "# x   Lambdax   Lambday   Lambdaz\n")
    # Write data in scientific format with 15 digits
    @inbounds begin
        for i = iL:iR
            write(f, @sprintf("%.15e", x[i]), "   ", @sprintf("%.15e", Λx[i]), "   ", @sprintf("%.15e", Λy[i]), "   ", @sprintf("%.15e", Λz[i]), "\n")
        end
    end
    # Close file
    close(f)
end

# Function to write out Reynolds stresses to a space delimited text file
function writeReynoldsStresses(t::AbstractFloat, x::AbstractArray, R11::AbstractArray, R22::AbstractArray, R33::AbstractArray, grid::rectilinearGrid, dataDir::AbstractString)
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
    @inbounds begin
        for i = iL:iR
            write(f, @sprintf("%.15e", x[i]), "   ", @sprintf("%.15e", R11[i]), "   ", @sprintf("%.15e", R22[i]), "   ", @sprintf("%.15e", R33[i]), "\n")
        end
    end
    # Close file
    close(f)
end

# Function to write out dissipation rates to a space delimited text file
function writeDissipationRates(t::AbstractFloat, x::AbstractArray, εx::AbstractArray, εy::AbstractArray, εz::AbstractArray, grid::rectilinearGrid, dataDir::AbstractString)
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
    @inbounds begin
        for i = iL:iR
            write(f, @sprintf("%.15e", x[i]), "   ", @sprintf("%.15e", εx[i]), "   ", @sprintf("%.15e", εy[i]), "   ", @sprintf("%.15e", εz[i]), "\n")
        end
    end
    # Close file
    close(f)
end

# Function to write out Taylor microscales to a space delimited text file
function writeTaylorMicroscales(t::AbstractFloat, x::AbstractArray, λx::AbstractArray, λy::AbstractArray, λz::AbstractArray, grid::rectilinearGrid, dataDir::AbstractString)
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
    @inbounds begin
        for i = iL:iR
            write(f, @sprintf("%.15e", x[i]), "   ", @sprintf("%.15e", λx[i]), "   ", @sprintf("%.15e", λy[i]), "   ", @sprintf("%.15e", λz[i]), "\n")
        end
    end
    # Close file
    close(f)
end

# Function to write out Kolmogorov microscales to a space delimited text file
function writeKolmogorovMicroscales(t::AbstractFloat, x::AbstractArray, ηx::AbstractArray, ηy::AbstractArray, ηz::AbstractArray, grid::rectilinearGrid, dataDir::AbstractString)
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
    @inbounds begin
        for i = iL:iR
            write(f, @sprintf("%.15e", x[i]), "   ", @sprintf("%.15e", ηx[i]), "   ", @sprintf("%.15e", ηy[i]), "   ", @sprintf("%.15e", ηz[i]), "\n")
        end
    end
    # Close file
    close(f)
end

# Function to write energy spectra to a space delimited text file
function writeEnergySpectra(t::AbstractFloat, κ::AbstractArray, Ex::AbstractArray, Ey::AbstractArray, Ez::AbstractArray, x::AbstractArray, grid::rectilinearGrid, dataDir::AbstractString)
    # Calculate Nyquist wavenumber
    N = Int(max(grid.Ny / 2, grid.Nz / 2))
    # Find index where x = xL
    iL = searchsortedfirst(x, grid.xL)
    # Find index where x = xR
    iR = searchsortedfirst(x, grid.xR)
    # Make file name for Ex
    filename = "$(dataDir)/spectraX_$(rpad(string(round(t, digits=5)), 7, "0")).dat"
    report("Writing energy spectra Ex to file $filename")
    # Open file
    f1 = open(filename, "w")
    # Write header
    write(f1, "# x   Ex\n")
    # Write data in scientific format with 15 digits
    @inbounds begin
        for i = iL:iR
            write(f1, @sprintf("%.15e", x[i]), "   ")
            for j = 1:N+1
                write(f1, @sprintf("%.15e", Ex[j, i]), "   ")
            end
            write(f1, "\n")
        end
    end
    # Close file
    close(f1)
    # Make file name for Ey
    filename = "$(dataDir)/spectraY_$(rpad(string(round(t, digits=5)), 7, "0")).dat"
    report("Writing energy spectra Ey to file $filename")
    # Open file
    f2 = open(filename, "w")
    # Write header
    write(f2, "# x   Ey\n")
    # Write data in scientific format with 15 digits
    @inbounds begin
        for i = iL:iR
            write(f2, @sprintf("%.15e", x[i]), "   ")
            for j = 1:N+1
                write(f2, @sprintf("%.15e", Ey[j, i]), "   ")
            end
            write(f2, "\n")
        end
    end
    # Close file
    close(f2)
    # Make file name for Ez
    filename = "$(dataDir)/spectraZ_$(rpad(string(round(t, digits=5)), 7, "0")).dat"
    report("Writing energy spectra Ez to file $filename")
    # Open file
    f3 = open(filename, "w")
    # Write header
    write(f3, "# x   Ez\n")
    # Write data in scientific format with 15 digits
    @inbounds begin
        for i = iL:iR
            write(f3, @sprintf("%.15e", x[i]), "   ")
            for j = 1:N+1
                write(f3, @sprintf("%.15e", Ez[j, i]), "   ")
            end
            write(f3, "\n")
        end
    end
    # Close file
    close(f3)
end

# Function to write integral length to a space delimited text file
function writeIntegralLength(t::AbstractFloat, Lyz::AbstractArray, x::AbstractArray, dataDir::AbstractString)
    # Find index where x = xL
    iL = searchsortedfirst(x, grid.xL)
    # Find index where x = xR
    iR = searchsortedfirst(x, grid.xR)
    # Get filename
    filename = "$(dataDir)/integralLength_$(rpad(string(round(t, digits=5)), 7, "0")).dat"
    report("Writing integral length to file $filename")
    # Open file
    f = open(filename, "w")
    # Write header
    write(f, "# x   Lyz\n")
    # Write data in scientific format with 15 digits
    @inbounds begin
        for i = iL:iR
            write(f, @sprintf("%.15e", x[i]), "   ", @sprintf("%.15e", Lyz[i]), "\n")
        end
    end
    # Close file
    close(f)
end