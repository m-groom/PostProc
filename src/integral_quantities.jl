# Functions for calculating and manipulating integral quantities

# Top level function to calculate integral quantities
function calcIntegralQuantities(t::Float64, x::SubArray{Float32,1}, y::SubArray{Float32,1}, z::SubArray{Float32,1}, Q::Array{Float32,4}, QBar::planeAverage, grid::rectilinearGrid, dataDir::String)
    # Calculate velocity correlation tensor
    tStart = report("Calculating velocity correlation tensor", 1)
    @time R11, R22, R33 = calcVelocityCorrelation(x, y, z, Q, QBar, grid, t, dataDir)
    tEnd = report("Finished calculating velocity correlation tensor...", 1)
    report("Elapsed time: $(tEnd - tStart)")
    exit()
    # Calculate correlation lengths
    tStart = report("Calculating correlation lengths", 1)
    calcCorrelationLengths(R11, R22, R33, x, y, z, grid)
    tEnd = report("Finished calculating correlation lengths...", 1)
    report("Elapsed time: $(tEnd - tStart)")
    # Calculate directional length scales
    tStart = report("Calculating directional length scales", 1)
    calcLengthScales(t, x, y, z, Q, QBar, grid, dataDir)
    tEnd = report("Finished calculating directional length scales...", 1)
    report("Elapsed time: $(tEnd - tStart)")
end

# Function to calculate velocity correlation tensor at x = x0
function calcVelocityCorrelation(x::SubArray{Float32,1}, y::SubArray{Float32,1}, z::SubArray{Float32,1}, Q::Array{Float32,4}, QBar::planeAverage, grid::rectilinearGrid, t::Float64, dataDir::String)
    # Find index where x = xL
    iL = searchsortedfirst(x, grid.xL)
    # Find index where x = xR
    iR = searchsortedfirst(x, grid.xR)
    # Get separation directions
    rx = @turbo x[iL:iR] .- 0.5 .* x[iR]
    ry = @turbo y .- 0.5 .* y[end]
    rz = @turbo z .- 0.5 .* z[end]
    # Initialise arrays
    R11 = zeros(Float64, length(rx), grid.Nx) # <u'(x)u'(x+rx)>
    R22 = zeros(Float64, length(ry), grid.Nx) # <v'(x)v'(x+ry)>
    R33 = zeros(Float64, length(rz), grid.Nx) # <w'(x)w'(x+rz)>
    nPtsInv = 1.0 / (grid.Ny * grid.Nz)
    # Loop over all cells
    @inbounds begin
        @batch for i = iL:iR
            for k = 1:grid.Nz
                for j = 1:grid.Ny
                    # Calculate fluctuating velocities at (i,j,k)
                    Ua = Q[i, j, k, 1] - QBar.UBar[i]
                    Va = Q[i, j, k, 2] - QBar.VBar[i]
                    Wa = Q[i, j, k, 3] - QBar.WBar[i]
                    # Calculate correlations in the x direction
                    for ii in eachindex(rx)
                        # Get index of point to calculate flucuation at
                        xp = x[i] + rx[ii]
                        # Interpolate velocity at xp
                        Ub = interp1D(xp, x, @view(Q[:, j, k, 1]), QBar.UBar)
                        # Calculate velocity correlation tensor components in the x direction
                        R11[ii, i] += Ua * Ub
                    end
                    # Calculate correlations in the y direction
                    for jj in eachindex(ry)
                        # Get index of point to calculate flucuation at
                        idx = Int(j + jj - grid.Ny / 2)
                        # Restrict y to be within 0 and 2π
                        if (idx <= 0)
                            idx += grid.Ny
                        elseif (idx > grid.Ny)
                            idx -= grid.Ny
                        end
                        # Interpolate velocity at idx
                        Vb = Q[i, idx, k, 2] - QBar.VBar[i]
                        # Calculate velocity correlation tensor components in the y direction
                        R22[jj, i] += Va * Vb
                    end
                    # Calculate correlations in the z direction
                    for kk in eachindex(rz)
                        # Get index of point to calculate flucuation at
                        idx = Int(k + kk - grid.Nz / 2)
                        # Restrict z to be within 0 and 2π
                        if (idx <= 0)
                            idx += grid.Nz
                        elseif (idx > grid.Nz)
                            idx -= grid.Nz
                        end
                        # Interpolate velocity at idx
                        Wb = Q[i, j, idx, 3] - QBar.WBar[i]
                        # Calculate velocity correlation tensor components in the z direction
                        R33[kk, i] += Wa * Wb
                    end
                end
            end
        end
    end
    # Divide by number of points to get averages
    @inbounds begin
        @tturbo for i = 1:grid.Nx
            for ii = 1:length(rx)
                R11[ii, i] *= nPtsInv
            end
        end
        @tturbo for i = 1:grid.Nx
            for jj = 1:length(ry)
                R22[jj, i] *= nPtsInv
            end
        end
        @tturbo for i = 1:grid.Nx
            for kk = 1:length(rz)
                R33[kk, i] *= nPtsInv
            end
        end
    end
    # Normalise by Rab(0)
    i0 = Int(ceil(length(rx) / 2))
    j0 = Int(ceil(length(ry) / 2))
    k0 = Int(ceil(length(rz) / 2))
    @inbounds begin
        @tturbo for i = iL:iR
            for ii = 1:length(rx)
                R11[ii, i] /= R11[i0, i]
            end
        end
        @tturbo for i = iL:iR
            for jj = 1:length(ry)
                R22[jj, i] /= R22[j0, i]
            end
        end
        @tturbo for i = iL:iR
            for kk = 1:length(rz)
                R33[kk, i] /= R33[k0, i]
            end
        end
    end
    # Write velocity correlations to file
    writeVelocityCorrelation(t, x, R11, R22, R33, grid, dataDir)

    return R11, R22, R33

end

# Function to calculate (longitudinal) correlation lengths
function calcCorrelationLengths(R11::Array{Float64,2}, R22::Array{Float64,2}, R33::Array{Float64,2}, x::SubArray{Float32,1}, y::SubArray{Float32,1}, z::SubArray{Float32,1}, grid::rectilinearGrid)
    # Find index where x = xL
    iL = searchsortedfirst(x, grid.xL)
    # Find index where x = xR
    iR = searchsortedfirst(x, grid.xR)
    # Get separation directions
    rx = @view(x[iL:iR]) .- 0.5 .* x[iR]
    ry = y .- 0.5 .* y[end]
    rz = z .- 0.5 .* z[end]
    # Get location of maximum
    i0 = Int(ceil(length(rx) / 2))
    j0 = Int(ceil(length(ry) / 2))
    k0 = Int(ceil(length(rz) / 2))
    # Initialise arrays
    Λx = zeros(Float64, grid.Nx)
    Λy = zeros(Float64, grid.Nx)
    Λz = zeros(Float64, grid.Nx)
    # Loop over cells in x direction
    @inbounds begin
        @simd for i = iL:iR
            # Get location of first zero crossing
            i2 = try
                i0 - 1 + findfirst(x -> x < 0.0, @view(R11[i0:end, i])) - 1
            catch
                size(R11, 1)
            end
            i1 = try
                findlast(x -> x < 0.0, @view(R11[1:i0, i])) + 1
            catch
                1
            end
            j2 = try
                j0 - 1 + findfirst(x -> x < 0.0, @view(R22[j0:end, i])) - 1
            catch
                size(R22, 1)
            end
            j1 = try
                findlast(x -> x < 0.0, @view(R22[1:j0, i])) + 1
            catch
                1
            end
            k2 = try
                k0 - 1 + findfirst(x -> x < 0.0, @view(R33[k0:end, i])) - 1
            catch
                size(R33, 1)
            end
            k1 = try
                findlast(x -> x < 0.0, @view(R33[1:k0, i])) + 1
            catch
                1
            end
            # Calculate correlation lengths
            Λx[i] = trapz(@view(rx[i1:i2]), @view(R11[i1:i2, i]))
            Λy[i] = trapz(@view(ry[j1:j2]), @view(R22[j1:j2, i]))
            Λz[i] = trapz(@view(rz[k1:k2]), @view(R33[k1:k2, i]))
        end
    end
    # Write correlation lengths to file
    writeCorrelationLengths(t, x, Λx, Λy, Λz, grid, dataDir)

end

# Function to calculate Taylor and Kolmorogov microscales
function calcLengthScales(t::Float64, x::SubArray{Float32,1}, y::SubArray{Float32,1}, z::SubArray{Float32,1}, Q::Array{Float32,4}, QBar::planeAverage, grid::rectilinearGrid, dataDir::String)
    # Initialise arrays
    R11 = zeros(Float64, grid.Nx) # <u'u'>
    R22 = zeros(Float64, grid.Nx) # <v'v'>
    R33 = zeros(Float64, grid.Nx) # <w'w'>
    dudxSquared = zeros(Float64, grid.Nx) # <(du/dx)^2>
    dvdySquared = zeros(Float64, grid.Nx) # <(dv/dy)^2>
    dwdzSquared = zeros(Float64, grid.Nx) # <(dw/dz)^2>
    omegaSquared = zeros(Float64, grid.Nx, 3) # <(omega_i)^2>
    divUSquared = zeros(Float64, grid.Nx) # <(du/dx+dv/dy+dw/dz)^2>
    nPtsInv = 1.0 / (grid.Ny * grid.Nz)
    # Loop over all cells
    @inbounds begin
        @batch for k = 1:grid.Nz
            for j = 1:grid.Ny
                @simd for i = 1:grid.Nx
                    uPrime = Q[i, j, k, 1] - QBar.UBar[i]
                    vPrime = Q[i, j, k, 2] - QBar.VBar[i]
                    wPrime = Q[i, j, k, 3] - QBar.WBar[i]
                    # Calculate Reynolds stresses
                    R11[i] += uPrime * uPrime
                    R22[i] += vPrime * vPrime
                    R33[i] += wPrime * wPrime
                    # Compute velocity derivatives
                    Uii = dUdX(x, @view(Q[:, j, k, 1]), i) # du/dx
                    Uji = dUdX(x, @view(Q[:, j, k, 2]), i) # dv/dx
                    Uki = dUdX(x, @view(Q[:, j, k, 3]), i) # dw/dx
                    Uij = dUdY(y, @view(Q[i, :, k, 1]), j) # du/dy
                    Ujj = dUdY(y, @view(Q[i, :, k, 2]), j) # dv/dy
                    Ukj = dUdY(y, @view(Q[i, :, k, 3]), j) # dw/dy
                    Uik = dUdY(z, @view(Q[i, j, :, 1]), k) # du/dz
                    Ujk = dUdY(z, @view(Q[i, j, :, 2]), k) # dv/dz
                    Ukk = dUdY(z, @view(Q[i, j, :, 3]), k) # dw/dz
                    # Compute square of velocity derivatives
                    dudxSquared[i] += Uii * Uii
                    dvdySquared[i] += Ujj * Ujj
                    dwdzSquared[i] += Ukk * Ukk
                    # Compute fluctuating vorticity
                    ωx = Ujk - Ukj
                    ωy = Uki - Uik
                    ωz = Uij - Uji
                    # Compute square of fluctuating vorticity
                    omegaSquared[i, 1] += ωx * ωx
                    omegaSquared[i, 2] += ωy * ωy
                    omegaSquared[i, 3] += ωz * ωz
                    # Compute fluctuating divergence
                    divU = Uii + Ujj + Ukk
                    # Compute square of fluctuating divergence
                    divUSquared = divU * divU
                end
            end
        end
    end
    # Divide by number of points to get averages
    @inbounds begin
        @tturbo for i = 1:grid.Nx
            R11[i] *= nPtsInv
            R22[i] *= nPtsInv
            R33[i] *= nPtsInv
            dudxSquared[i] *= nPtsInv
            dvdySquared[i] *= nPtsInv
            dwdzSquared[i] *= nPtsInv
            omegaSquared[i] *= nPtsInv
            divUSquared[i] *= nPtsInv
        end
    end
    # Calculate dissipation rates
    εx, εy, εz = calcDissipationRates(QBar, omegaSquared, divUSquared, grid.Nx)
    # Calculate Taylor microscales
    λx, λy, λz = calcTaylorMicroscales(R11, R22, R33, dudxSquared, dvdySquared, dwdzSquared, grid.Nx)
    # Calculate Kolmogorov microscales
    ηx, ηy, ηz = calcKolmogorovMicroscales(QBar.nuBar, εx, εy, εz, grid.Nx)
    # Write Reynolds stresses to file
    writeReynoldsStresses(t, x, R11, R22, R33, grid, dataDir)
    # Write dissipation rates to file
    writeDissipationRates(t, x, εx, εy, εz, grid, dataDir)
    # Write Taylor microscales to file
    writeTaylorMicroscales(t, x, λx, λy, λz, grid, dataDir)
    # Write Kolmogorov microscales to file
    writeKolmogorovMicroscales(t, x, ηx, ηy, ηz, grid, dataDir)

end

# Function to calculate dissipation rates
function calcDissipationRates(QBar::planeAverage, omegaSquared::Array{Float64,2}, divUSquared::Array{Float64,1}, Nx::Int64)
    # Initialise arrays
    εx = zeros(Float64, Nx)
    εy = zeros(Float64, Nx)
    εz = zeros(Float64, Nx)
    # Store 4/9 as a constant
    fourNinths = 4.0 / 9.0
    # Loop over all cells
    @inbounds begin
        @tturbo for i = 1:Nx
            nuBar = QBar.muBar[i] / QBar.rhoBar[i]
            εx[i] = nuBar * (omegaSquared[i, 1] + fourNinths * divUSquared[i])
            εy[i] = nuBar * (omegaSquared[i, 2] + fourNinths * divUSquared[i])
            εz[i] = nuBar * (omegaSquared[i, 3] + fourNinths * divUSquared[i])
        end
    end
    return εx, εy, εz
end

# Function to calculate Taylor microscales
function calcTaylorMicroscales(R11::Array{Float64,1}, R22::Array{Float64,1}, R33::Array{Float64,1}, dudxSquared::Array{Float64,1}, dvdySquared::Array{Float64,1}, dwdzSquared::Array{Float64,1}, Nx::Int64)
    # Initialise arrays
    λx = zeros(Float64, Nx)
    λy = zeros(Float64, Nx)
    λz = zeros(Float64, Nx)
    # Loop over all cells
    @inbounds begin
        @tturbo for i = 1:Nx
            λx[i] = sqrt(R11[i] / dudxSquared[i])
            λy[i] = sqrt(R22[i] / dvdySquared[i])
            λz[i] = sqrt(R33[i] / dwdzSquared[i])
        end
    end
    return λx, λy, λz
end

# Function to calculate Kolmogorov microscales
function calcKolmogorovMicroscales(nuBar::Array{Float64,1}, εx::Array{Float64,1}, εy::Array{Float64,1}, εz::Array{Float64,1}, Nx::Int64)
    # Initialise arrays
    ηx = zeros(Float64, Nx)
    ηy = zeros(Float64, Nx)
    ηz = zeros(Float64, Nx)
    # Loop over all cells
    @inbounds begin
        @tturbo for i = 1:Nx
            nuBarCubed = nuBar[i]^3
            ηx[i] = (nuBarCubed / εx[i])^0.25
            ηy[i] = (nuBarCubed / εy[i])^0.25
            ηz[i] = (nuBarCubed / εz[i])^0.25
        end
    end
    return ηx, ηy, ηz
end