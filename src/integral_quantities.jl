# Functions for calculating and manipulating integral quantities

# Top level function to calculate integral quantities
function calcIntegralQuantities(t::Float64, x::SubArray{Float32,1}, y::SubArray{Float32,1}, z::SubArray{Float32,1}, Q::Array{Float32,4}, QBar::planeAverage, grid::rectilinearGrid, nVars::Int64, x0::Float64, dataDir::String)
    # Calculate velocity correlation tensor
    tStart = report("Calculating velocity correlation tensor", 1)
    R11, R22, R33 = calcVelocityCorrelation(x, y, z, Q, QBar, grid, nVars, x0)
    tEnd = report("Finished calculating velocity correlation tensor...", 1)
    report("Elapsed time: $(tEnd - tStart)")
    # Calculate correlation lengths
    tStart = report("Calculating correlation lengths", 1)
    Λx, Λyz = calcCorrelationLengths(R11, R22, R33, x, y, z, grid)
    tEnd = report("Finished calculating correlation lengths...", 1)
    report("Elapsed time: $(tEnd - tStart)")
    # Write correlation lengths to file
    writeCorrelationLengths(t, Λx, Λyz, dataDir)
    # Calculate directional length scales
    tStart = report("Calculating directional length scales", 1)
    λx, λyz, ηx, ηyz = calcLengthScales(t, x, y, z, Q, QBar, grid, nVars, x0, dataDir)
    tEnd = report("Finished calculating directional length scales...", 1)
    report("Elapsed time: $(tEnd - tStart)")
    # Write directional length scales to file
    writeLengthScales(t, λx, λyz, ηx, ηyz, dataDir)
    # Calculate integral width
    tStart = report("Calculating integral width", 1)
    W, H = calcIntegralWidth(t, QBar, dataDir)
    tEnd = report("Finished calculating integral width...", 1)
    report("Elapsed time: $(tEnd - tStart)")
end

# Function to calculate velocity correlation tensor at x = x0
function calcVelocityCorrelation(x::SubArray{Float32,1}, y::SubArray{Float32,1}, z::SubArray{Float32,1}, Q::Array{Float32,4}, QBar::planeAverage, grid::rectilinearGrid, nVars::Int64, x0::Float64)
    # Find index where x = xL
    iL = searchsortedfirst(x, grid.xL)
    # Find index where x = xR
    iR = searchsortedfirst(x, grid.xR)
    # Get grid extents
    yextent = grid.yR - grid.yL
    zextent = grid.zR - grid.zL
    # Get separation directions
    rx = @view(x[iL:iR]) .- 0.5 .* x[iR]
    ry = y .- 0.5 .* y[end]
    rz = z .- 0.5 .* z[end]
    # Initialise arrays
    R11 = zeros(Float64, length(rx)) # <u'(x)u'(x+rx)>
    R22 = zeros(Float64, length(ry)) # <v'(x)v'(x+ry)>
    R33 = zeros(Float64, length(rz)) # <w'(x)w'(x+rz)>
    nPtsInv = 1.0 / (grid.Ny * grid.Nz)
    # Find index where x = x0
    i = searchsortedfirst(x, x0) - 1
    # Loop over all cells
    @inbounds for k = 1:grid.Nz
        @inbounds for j = 1:grid.Ny
            # Calculate fluctuating velocities at (i,j,k)
            DaInv = 1.0 / Q[i, j, k, nVars-3]
            Ua = Q[i, j, k, 1] * DaInv - QBar.UBar[i]
            Va = Q[i, j, k, 2] * DaInv - QBar.VBar[i]
            Wa = Q[i, j, k, 3] * DaInv - QBar.WBar[i]
            # Calculate correlations in the x direction
            @inbounds for ii in eachindex(rx)
                # Get index of point to calculate flucuation at
                xp = x[i] + rx[ii]
                idx = searchsortedfirst(x, xp) - 1
                # Calculate fluctuating velocities at (idx,j,k)
                DbInv = 1.0 / Q[idx, j, k, nVars-3]
                Ub = Q[idx, j, k, 1] * DbInv - QBar.UBar[idx]
                # Calculate velocity correlation tensor components in the x direction
                R11[ii] += Ua * Ub
            end
            # Calculate correlations in the y direction
            @inbounds for jj in eachindex(ry)
                # Get index of point to calculate flucuation at
                yp = y[j] + ry[jj]
                # Restrict y to be within 0 and 2π
                if yp < grid.yL
                    yp += yextent
                elseif yp > grid.yR
                    yp -= yextent
                end
                idx = searchsortedfirst(y, yp) - 1
                # Calculate fluctuating velocities at (i,idx,k)
                DbInv = 1.0 / Q[i, idx, k, nVars-3]
                Vb = Q[i, idx, k, 2] * DbInv - QBar.VBar[i]
                # Calculate velocity correlation tensor components in the y direction
                R22[jj] += Va * Vb
            end
            # Calculate correlations in the z direction
            @inbounds for kk in eachindex(rz)
                # Get index of point to calculate flucuation at
                zp = z[k] + rz[kk]
                # Restrict z to be within 0 and 2π
                if zp < grid.zL
                    zp += zextent
                elseif zp > grid.zR
                    zp -= zextent
                end
                idx = searchsortedfirst(z, zp) - 1
                # Calculate fluctuating velocities at (i,j,idx)
                DbInv = 1.0 / Q[i, j, idx, nVars-3]
                Wb = Q[i, j, idx, 3] * DbInv - QBar.WBar[i]
                # Calculate velocity correlation tensor components in the z direction
                R33[kk] += Wa * Wb
            end
        end
    end
    # Divide by number of points to get averages
    R11 *= nPtsInv
    R22 *= nPtsInv
    R33 *= nPtsInv
    # Normalise by Rab(0)
    i0 = Int(ceil(length(rx) / 2))
    R11 /= R11[i0]
    j0 = Int(ceil(length(ry) / 2))
    R22 /= R22[j0]
    k0 = Int(ceil(length(rz) / 2))
    R33 /= R33[k0]

    return R11, R22, R33

end

# Function to calculate (longitudinal) correlation lengths
function calcCorrelationLengths(R11::Array{Float64,1}, R22::Array{Float64,1}, R33::Array{Float64,1}, x::SubArray{Float32,1}, y::SubArray{Float32,1}, z::SubArray{Float32,1}, grid::rectilinearGrid)
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
    # Get location of first zero crossing
    i2 = try
        i0 - 1 + findfirst(x -> x < 0.0, @view(R11[i0:end])) - 1
    catch
        length(R11)
    end
    i1 = try
        findlast(x -> x < 0.0, @view(R11[1:i0])) + 1
    catch
        1
    end
    j2 = try
        j0 - 1 + findfirst(x -> x < 0.0, @view(R22[j0:end])) - 1
    catch
        length(R22)
    end
    j1 = try
        findlast(x -> x < 0.0, @view(R22[1:j0])) + 1
    catch
        1
    end
    k2 = try
        k0 - 1 + findfirst(x -> x < 0.0, @view(R33[k0:end])) - 1
    catch
        length(R33)
    end
    k1 = try
        findlast(x -> x < 0.0, @view(R33[1:k0])) + 1
    catch
        1
    end
    # Calculate correlation lengths
    Λx = trapz(@view(rx[i1:i2]), @view(R11[i1:i2]))
    Λy = trapz(@view(ry[j1:j2]), @view(R22[j1:j2]))
    Λz = trapz(@view(rz[k1:k2]), @view(R33[k1:k2]))
    Λyz = 0.5 .* (Λy .+ Λz)

    return Λx, Λyz

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
    λyz = 0.5 .* (λy .+ λz)
    ηx = (QBar.nuBar .^ 3 ./ εx) .^ 0.25
    ηy = (QBar.nuBar .^ 3 ./ εy) .^ 0.25
    ηz = (QBar.nuBar .^ 3 ./ εz) .^ 0.25
    ηyz = 0.5 .* (ηy .+ ηz)
    # Write Reynolds stresses to file
    writeReynoldsStresses(t, x, R11, R22, R33, grid, dataDir)
    # Write dissipation rates to file
    writeDissipationRates(t, x, εx, εy, εz, grid, dataDir)
    # Write Taylor microscales to file
    writeTaylorMicroscales(t, x, λx, λy, λz, grid, dataDir)
    # Write Kolmogorov microscales to file
    writeKolmogorovMicroscales(t, x, ηx, ηy, ηz, grid, dataDir)
    # Get location of interface
    i0 = searchsortedfirst(x, x0) - 1

    return λx[i0], λyz[i0], ηx[i0], ηyz[i0]

end