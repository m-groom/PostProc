# Functions for calculating and manipulating integral quantities

# Top level function to calculate integral quantities
function calcIntegralQuantities(t::Float64, x::SubArray{Float32,1}, y::SubArray{Float32,1}, z::SubArray{Float32,1}, Q::Array{Float32,4}, QBar::planeAverage, grid::rectilinearGrid, nVars::Int64, x0::Float64, dataDir::String)
    # Calculate velocity correlation tensor
    tStart = report("Calculating velocity correlation tensor", 1)
    R11, R22, R33 = calcVelocityCorrelation(x, y, z, Q, QBar, grid, t, dataDir)
    tEnd = report("Finished calculating velocity correlation tensor...", 1)
    report("Elapsed time: $(tEnd - tStart)")
    # Calculate correlation lengths
    tStart = report("Calculating correlation lengths", 1)
    calcCorrelationLengths(R11, R22, R33, x, y, z, grid)
    tEnd = report("Finished calculating correlation lengths...", 1)
    report("Elapsed time: $(tEnd - tStart)")
    # Calculate directional length scales
    tStart = report("Calculating directional length scales", 1)
    @time calcLengthScales(t, x, y, z, Q, QBar, grid, nVars, x0, dataDir)
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
    rx = @view(x[iL:iR]) .- 0.5 .* x[iR]
    ry = y .- 0.5 .* y[end]
    rz = z .- 0.5 .* z[end]
    # Initialise arrays
    R11 = zeros(Float64, length(rx), grid.Nx) # <u'(x)u'(x+rx)>
    R22 = zeros(Float64, length(ry), grid.Nx) # <v'(x)v'(x+ry)>
    R33 = zeros(Float64, length(rz), grid.Nx) # <w'(x)w'(x+rz)>
    nPtsInv = 1.0 / (grid.Ny * grid.Nz)
    # Loop over all cells
    @inbounds for k = 1:grid.Nz
        @inbounds for j = 1:grid.Ny
            @inbounds for i = iL:iR
                # Calculate fluctuating velocities at (i,j,k)
                Ua = Q[i, j, k, 1] - QBar.UBar[i]
                Va = Q[i, j, k, 2] - QBar.VBar[i]
                Wa = Q[i, j, k, 3] - QBar.WBar[i]
                # Calculate correlations in the x direction
                @inbounds for ii in eachindex(rx)
                    # Get index of point to calculate flucuation at
                    xp = x[i] + rx[ii]
                    # Interpolate velocity at xp
                    Ub = interp1D(xp, x, @view(Q[:, j, k, 1]), QBar.UBar)
                    # Calculate velocity correlation tensor components in the x direction
                    R11[ii, i] += Ua * Ub
                end
                # Calculate correlations in the y direction
                @inbounds for jj in eachindex(ry)
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
                @inbounds for kk in eachindex(rz)
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
    # Divide by number of points to get averages
    R11 *= nPtsInv
    R22 *= nPtsInv
    R33 *= nPtsInv
    # Normalise by Rab(0)
    i0 = Int(ceil(length(rx) / 2))
    j0 = Int(ceil(length(ry) / 2))
    k0 = Int(ceil(length(rz) / 2)) 
    @inbounds for i = iL:iR
        R11[:, i] /= R11[i0, i]
        R22[:, i] /= R22[j0, i]
        R33[:, i] /= R33[k0, i]
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
    @inbounds for i = iL:iR
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
    # Write correlation lengths to file
    writeCorrelationLengths(t, x, Λx, Λy, Λz, grid, dataDir)

end

# Function to calculate Taylor and Kolmorogov microscales
function calcLengthScales(t::Float64, x::SubArray{Float32,1}, y::SubArray{Float32,1}, z::SubArray{Float32,1}, Q::Array{Float32,4}, QBar::planeAverage, grid::rectilinearGrid, nVars::Int64, x0::Float64, dataDir::String)
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
    @inbounds for k = 1:grid.Nz
        @inbounds for j = 1:grid.Ny
            @inbounds for i = 1:grid.Nx
                rhoInv = 1.0 / Q[i, j, k, nVars-3]
                uPrime = Q[i, j, k, 1] * rhoInv - QBar.UBar[i]
                vPrime = Q[i, j, k, 2] * rhoInv - QBar.VBar[i]
                wPrime = Q[i, j, k, 3] * rhoInv - QBar.WBar[i]
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
    # Divide by number of points to get averages
    R11 *= nPtsInv
    R22 *= nPtsInv
    R33 *= nPtsInv
    dudxSquared *= nPtsInv
    dvdySquared *= nPtsInv
    dwdzSquared *= nPtsInv
    omegaSquared *= nPtsInv
    divUSquared *= nPtsInv
    # Calculate dissipation rate
    εx = QBar.muBar ./ QBar.rhoBar .* (omegaSquared[:, 1] .+ 4.0 / 9.0 .* divUSquared)
    εy = QBar.muBar ./ QBar.rhoBar .* (omegaSquared[:, 2] .+ 4.0 / 9.0 .* divUSquared)
    εz = QBar.muBar ./ QBar.rhoBar .* (omegaSquared[:, 3] .+ 4.0 / 9.0 .* divUSquared)
    # Calculate Taylor and Kolmogorov microscales
    λx = sqrt.(R11 ./ dudxSquared)
    λy = sqrt.(R22 ./ dvdySquared)
    λz = sqrt.(R33 ./ dwdzSquared)
    ηx = (QBar.nuBar .^ 3 ./ εx) .^ 0.25
    ηy = (QBar.nuBar .^ 3 ./ εy) .^ 0.25
    ηz = (QBar.nuBar .^ 3 ./ εz) .^ 0.25  
    # Write Reynolds stresses to file
    writeReynoldsStresses(t, x, R11, R22, R33, grid, dataDir)
    # Write dissipation rates to file
    writeDissipationRates(t, x, εx, εy, εz, grid, dataDir)
    # Write Taylor microscales to file
    writeTaylorMicroscales(t, x, λx, λy, λz, grid, dataDir)
    # Write Kolmogorov microscales to file
    writeKolmogorovMicroscales(t, x, ηx, ηy, ηz, grid, dataDir)

end