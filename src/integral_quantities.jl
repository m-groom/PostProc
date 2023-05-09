# Functions for calculating and manipulating integral quantities

# Top level function to calculate integral quantities
function calcIntegralQuantities(t::Float64, x::Array{Float32,3}, y::Array{Float32,3}, z::Array{Float32,3}, Q::Array{Float32,4}, QBar::planeAverage, grid::rectilinearGrid, nVars::Int64, x0::Float64)
    # Calculate velocity correlation tensor
    tStart = report("Calculating velocity correlation tensor", 1)
    R11, R22, R33 = calcVelocityCorrelation(x[:, 1, 1], y[1, :, 1], z[1, 1, :], Q, QBar, grid, nVars, x0)
    tEnd = report("Finished calculating velocity correlation tensor...", 1)
    report("Elapsed time: $(tEnd - tStart)")
    # Calculate correlation lengths
    tStart = report("Calculating correlation lengths", 1)
    Λx, Λyz = calcCorrelationLengths(R11, R22, R33, x[:, 1, 1], y[1, :, 1], z[1, 1, :], grid)
    tEnd = report("Finished calculating correlation lengths...", 1)
    report("Elapsed time: $(tEnd - tStart)")
    # Write correlation lengths to file
    writeCorrelationLengths(t, Λx, Λyz)
    # Calculate directional length scales
    tStart = report("Calculating directional length scales", 1)
    λx, λyz, ηx, ηyz = calcLengthScales(t, x[:, 1, 1], y[1, :, 1], z[1, 1, :], Q, QBar, grid, nVars, x0)
    tEnd = report("Finished calculating directional length scales...", 1)
    report("Elapsed time: $(tEnd - tStart)")
    # Write directional length scales to file
    writeLengthScales(t, λx, λyz, ηx, ηyz)
    # Calculate integral width
    tStart = report("Calculating integral width", 1)
    W, H, Wb, Ws, Hb, Hs = calcIntegralWidth(t, QBar, x0)
    tEnd = report("Finished calculating integral width...", 1)
    report("Elapsed time: $(tEnd - tStart)")
end

# Function to calculate velocity correlation tensor at x = x0
function calcVelocityCorrelation(x::Array{Float32,1}, y::Array{Float32,1}, z::Array{Float32,1}, Q::Array{Float32,4}, QBar::planeAverage, grid::rectilinearGrid, nVars::Int64, x0::Float64)
    # Find index where x = xL
    iL = argmin(abs.(x .- grid.xL))
    # Find index where x = xR
    iR = argmin(abs.(x .- grid.xR))
    # Get separation directions
    rx = x[iL:iR] .- 0.5 * x[iR]
    ry = y .- 0.5 * y[end]
    rz = z .- 0.5 * z[end]
    # Initialise arrays
    R11 = zeros(Float64, length(rx)) # <u'(x)u'(x+rx)>
    # R21 = zeros(Float64, length(rx)) # <v'(x)v'(x+rx)>
    # R31 = zeros(Float64, length(rx)) # <w'(x)w'(x+rx)>
    # R12 = zeros(Float64, length(ry)) # <u'(x)u'(x+ry)>
    R22 = zeros(Float64, length(ry)) # <v'(x)v'(x+ry)>
    # R32 = zeros(Float64, length(ry)) # <w'(x)w'(x+ry)>
    # R13 = zeros(Float64, length(rz)) # <u'(x)u'(x+rz)>
    # R23 = zeros(Float64, length(rz)) # <v'(x)v'(x+rz)>
    R33 = zeros(Float64, length(rz)) # <w'(x)w'(x+rz)>
    nPtsInv = 1.0 / (grid.Ny * grid.Nz)
    # Find index where x = x0
    i = argmin(abs.(x .- x0))
    # Loop over all cells
    for k = 1:grid.Nz
        for j = 1:grid.Ny
            # Calculate fluctuating velocities at (i,j,k)
            DaInv = 1.0 / Q[i, j, k, nVars-3]
            Ua = Q[i, j, k, 1] * DaInv - QBar.UBar[i]
            Va = Q[i, j, k, 2] * DaInv - QBar.VBar[i]
            Wa = Q[i, j, k, 3] * DaInv - QBar.WBar[i]
            # Calculate correlations in the x direction
            for ii in eachindex(rx)
                # Get index of point to calculate flucuation at
                xp = x[i] + rx[ii]
                idx = argmin(abs.(x .- xp))
                # Calculate fluctuating velocities at (idx,j,k)
                DbInv = 1.0 / Q[idx, j, k, nVars-3]
                Ub = Q[idx, j, k, 1] * DbInv - QBar.UBar[idx]
                # Vb = Q[idx, j, k, 2] * DbInv - QBar.VBar[idx]
                # Wb = Q[idx, j, k, 3] * DbInv - QBar.WBar[idx]
                # Calculate velocity correlation tensor components in the x direction
                R11[ii] = R11[ii] + Ua * Ub
                # R21[ii] = R21[ii] + Va * Vb
                # R31[ii] = R31[ii] + Wa * Wb
            end
            # Calculate correlations in the y direction
            for jj in eachindex(ry)
                # Get index of point to calculate flucuation at
                yp = y[j] + ry[jj]
                # Restrict y to be within 0 and 2π
                if yp < grid.yL
                    yp = yp + (grid.yR - grid.yL)
                elseif yp > grid.yR
                    yp = yp - (grid.yR - grid.yL)
                end
                idx = argmin(abs.(y .- yp))
                # Calculate fluctuating velocities at (i,idx,k)
                DbInv = 1.0 / Q[i, idx, k, nVars-3]
                # Ub = Q[i, idx, k, 1] * DbInv - QBar.UBar[i]
                Vb = Q[i, idx, k, 2] * DbInv - QBar.VBar[i]
                # Wb = Q[i, idx, k, 3] * DbInv - QBar.WBar[i]
                # Calculate velocity correlation tensor components in the y direction
                # R12[jj] = R12[jj] + Ua * Ub
                R22[jj] = R22[jj] + Va * Vb
                # R32[jj] = R32[jj] + Wa * Wb
            end
            # Calculate correlations in the z direction
            for kk in eachindex(rz)
                # Get index of point to calculate flucuation at
                zp = z[k] + rz[kk]
                # Restrict z to be within 0 and 2π
                if zp < grid.zL
                    zp = zp + (grid.zR - grid.zL)
                elseif zp > grid.zR
                    zp = zp - (grid.zR - grid.zL)
                end
                idx = argmin(abs.(z .- zp))
                # Calculate fluctuating velocities at (i,j,idx)
                DbInv = 1.0 / Q[i, j, idx, nVars-3]
                # Ub = Q[i, j, idx, 1] * DbInv - QBar.UBar[i]
                # Vb = Q[i, j, idx, 2] * DbInv - QBar.VBar[i]
                Wb = Q[i, j, idx, 3] * DbInv - QBar.WBar[i]
                # Calculate velocity correlation tensor components in the z direction
                # R13[kk] = R13[kk] + Ua * Ub
                # R23[kk] = R23[kk] + Va * Vb
                R33[kk] = R33[kk] + Wa * Wb
            end
        end
    end
    # Divide by number of points to get averages
    R11 = R11 * nPtsInv
    # R21 = R21 * nPtsInv
    # R31 = R31 * nPtsInv
    # R12 = R12 * nPtsInv
    R22 = R22 * nPtsInv
    # R32 = R32 * nPtsInv
    # R13 = R13 * nPtsInv
    # R23 = R23 * nPtsInv
    R33 = R33 * nPtsInv
    # Normalise by Rab(0)
    i0 = Int(ceil(length(rx) / 2))
    R11 = R11 / R11[i0]
    # R21 = R21 / R21[i0]
    # R31 = R31 / R31[i0]
    j0 = Int(ceil(length(ry) / 2))
    # R12 = R12 / R12[j0]
    R22 = R22 / R22[j0]
    # R32 = R32 / R32[j0]
    k0 = Int(ceil(length(rz) / 2))
    # R13 = R13 / R13[k0]
    # R23 = R23 / R23[k0]
    R33 = R33 / R33[k0]

    return R11, R22, R33

end

# Function to calculate (longitudinal) correlation lengths
function calcCorrelationLengths(R11::Array{Float64,1}, R22::Array{Float64,1}, R33::Array{Float64,1}, x::Array{Float32,1}, y::Array{Float32,1}, z::Array{Float32,1}, grid::rectilinearGrid)
    # Find index where x = xL
    iL = argmin(abs.(x .- grid.xL))
    # Find index where x = xR
    iR = argmin(abs.(x .- grid.xR))
    # Get separation directions
    rx = x[iL:iR] .- 0.5 * x[iR]
    ry = y .- 0.5 * y[end]
    rz = z .- 0.5 * z[end]
    # Get location of maximum
    i0 = Int(ceil(length(rx) / 2))
    j0 = Int(ceil(length(ry) / 2))
    k0 = Int(ceil(length(rz) / 2))
    # Get location of first zero crossing
    if (all(x -> x >= 0.0, R11)) # Check if R11 is positive everywhere
        i1 = 1
        i2 = length(R11)
    else
        i2 = i0 - 1 + findfirst(x -> x < 0.0, R11[i0:end]) - 1
        i1 = findlast(x -> x < 0.0, R11[1:i0]) + 1
    end
    if (all(x -> x >= 0.0, R22)) # Check if R22 is positive everywhere
        j1 = 1
        j2 = length(R22)
    else
        j2 = j0 - 1 + findfirst(x -> x < 0.0, R22[j0:end]) - 1
        j1 = findlast(x -> x < 0.0, R22[1:j0]) + 1
    end
    if (all(x -> x >= 0.0, R33)) # Check if R33 is positive everywhere
        k1 = 1
        k2 = length(R33)
    else
        k2 = k0 - 1 + findfirst(x -> x < 0.0, R33[k0:end]) - 1
        k1 = findlast(x -> x < 0.0, R33[1:k0]) + 1
    end
    # Calculate correlation lengths
    Λx = trapz(rx[i1:i2], R11[i1:i2])
    Λy = trapz(ry[j1:j2], R22[j1:j2])
    Λz = trapz(rz[k1:k2], R33[k1:k2])
    Λyz = 0.5 * (Λy + Λz)

    return Λx, Λyz

end

# Function to calculate Taylor and Kolmorogov microscales
function calcLengthScales(t::Float64, x::Array{Float32,1}, y::Array{Float32,1}, z::Array{Float32,1}, Q::Array{Float32,4}, QBar::planeAverage, grid::rectilinearGrid, nVars::Int64, x0::Float64)
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
    for k = 1:grid.Nz
        for j = 1:grid.Ny
            for i = 1:grid.Nx
                rhoInv = 1.0 / Q[i, j, k, nVars-3]
                u = Q[i, j, k, 1] * rhoInv - QBar.UBar[i]
                v = Q[i, j, k, 2] * rhoInv - QBar.VBar[i]
                w = Q[i, j, k, 3] * rhoInv - QBar.WBar[i]
                # Calculate Reynolds stresses
                R11[i] = R11[i] + u * u
                R22[i] = R22[i] + v * v
                R33[i] = R33[i] + w * w
                # Compute velocity derivatives
                Uii = dUdX(x, Q[:, j, k, 1], i) # du/dx
                Uji = dUdX(x, Q[:, j, k, 2], j) # dv/dx
                Uki = dUdX(x, Q[:, j, k, 3], i) # dw/dx
                Uij = dUdY(y, Q[i, :, k, 1], j) # du/dy
                Ujj = dUdY(y, Q[i, :, k, 2], j) # dv/dy
                Ukj = dUdY(y, Q[i, :, k, 3], j) # dw/dy
                Uik = dUdY(z, Q[i, j, :, 1], k) # du/dz
                Ujk = dUdY(z, Q[i, j, :, 2], k) # dv/dz
                Ukk = dUdY(z, Q[i, j, :, 3], k) # dw/dz
                # Compute square of velocity derivatives
                dudxSquared[i] = dudxSquared[i] + Uii * Uii
                dvdySquared[i] = dvdySquared[i] + Ujj * Ujj
                dwdzSquared[i] = dwdzSquared[i] + Ukk * Ukk
                # Compute fluctuating vorticity
                ωx = Ujk - Ukj
                ωy = Uki - Uik
                ωz = Uij - Uji
                # Compute square of fluctuating vorticity
                omegaSquared[i, 1] = omegaSquared[i, 1] + ωx * ωx
                omegaSquared[i, 2] = omegaSquared[i, 2] + ωy * ωy
                omegaSquared[i, 3] = omegaSquared[i, 3] + ωz * ωz
                # Compute fluctuating divergence
                divU = Uii + Ujj + Ukk
                # Compute square of fluctuating divergence
                divUSquared = divU * divU
            end
        end
    end
    # Divide by number of points to get averages
    R11 = R11 * nPtsInv
    R22 = R22 * nPtsInv
    R33 = R33 * nPtsInv
    dudxSquared = dudxSquared * nPtsInv
    dvdySquared = dvdySquared * nPtsInv
    dwdzSquared = dwdzSquared * nPtsInv
    omegaSquared = omegaSquared * nPtsInv
    divUSquared = divUSquared * nPtsInv
    # Calculate dissipation rate
    εx = QBar.muBar ./ QBar.rhoBar .* (omegaSquared[:, 1] .+ 4.0 / 9.0 * divUSquared)
    εy = QBar.muBar ./ QBar.rhoBar .* (omegaSquared[:, 2] .+ 4.0 / 9.0 * divUSquared)
    εz = QBar.muBar ./ QBar.rhoBar .* (omegaSquared[:, 3] .+ 4.0 / 9.0 * divUSquared)
    # Calculate Taylor and Kolmogorov microscales
    λx = sqrt.(R11 ./ dudxSquared)
    λy = sqrt.(R22 ./ dvdySquared)
    λz = sqrt.(R33 ./ dwdzSquared)
    λyz = 0.5 * (λy + λz)
    ηx = (QBar.nuBar .^ 3 ./ εx) .^ 0.25
    ηy = (QBar.nuBar .^ 3 ./ εy) .^ 0.25
    ηz = (QBar.nuBar .^ 3 ./ εz) .^ 0.25
    ηyz = 0.5 * (ηy + ηz)
    # Write Reynolds stresses to file
    writeReynoldsStresses(t, x, R11, R22, R33)
    # Write dissipation rates to file
    writeDissipationRates(t, x, εx, εy, εz)
    # Write Taylor microscales to file
    writeTaylorMicroscales(t, x, λx, λy, λz)
    # Write Kolmogorov microscales to file
    writeKolmogorovMicroscales(t, x, ηx, ηy, ηz)
    # Get location of interface
    i0 = argmin(abs.(x .- x0))

    return λx[i0], λyz[i0], ηx[i0], ηyz[i0]

end