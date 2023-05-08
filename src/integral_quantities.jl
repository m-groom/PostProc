# Functions for calculating and manipulating integral quantities

# Top level function to calculate integral quantities
function calcIntegralQuantities(t::Float64, x::Array{Float32,3}, y::Array{Float32,3}, z::Array{Float32,3}, Q::Array{Float32,4}, QBar::planeAverage, grid::rectilinearGrid, nVars::Int64, x0::Float64)
    # Calculate velocity correlation tensor
    tStart = report("Calculating velocity correlation tensor", 1)
    R11, R21, R31, R12, R22, R32, R13, R23, R33 = calcVelocityCorrelation(x[:, 1, 1], y[1, :, 1], z[1, 1, :], Q, QBar, grid, nVars, x0)
    tEnd = report("Finished calculating velocity correlation tensor...", 1)
    report("Elapsed time: $(tEnd - tStart)")
    # Calculate correlation lengths
    Λx, Λyz = calcCorrelationLengths(R11, R22, R33, x[:, 1, 1], y[1, :, 1], z[1, 1, :], grid)

    # Calculate directional length scales
    # λx, λyz, ηx, ηyz = calcLengthScales(x, y, z, Q, QBar, grid, nVars)

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
    R11 = zeros(Float64, length(rx))
    R21 = zeros(Float64, length(rx))
    R31 = zeros(Float64, length(rx))
    R12 = zeros(Float64, length(ry))
    R22 = zeros(Float64, length(ry))
    R32 = zeros(Float64, length(ry))
    R13 = zeros(Float64, length(rz))
    R23 = zeros(Float64, length(rz))
    R33 = zeros(Float64, length(rz))
    nPtsInv = 1.0 / (grid.Ny * grid.Nz)
    # Find index where x = x0
    i = argmin(abs.(x[:, 1, 1] .- x0))
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
                Vb = Q[idx, j, k, 2] * DbInv - QBar.VBar[idx]
                Wb = Q[idx, j, k, 3] * DbInv - QBar.WBar[idx]
                # Calculate velocity correlation tensor components in the x direction
                R11[ii] = R11[ii] + Ua * Ub
                R21[ii] = R21[ii] + Va * Vb
                R31[ii] = R31[ii] + Wa * Wb
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
                Ub = Q[i, idx, k, 1] * DbInv - QBar.UBar[i]
                Vb = Q[i, idx, k, 2] * DbInv - QBar.VBar[i]
                Wb = Q[i, idx, k, 3] * DbInv - QBar.WBar[i]
                # Calculate velocity correlation tensor components in the y direction
                R12[jj] = R12[jj] + Ua * Ub
                R22[jj] = R22[jj] + Va * Vb
                R32[jj] = R32[jj] + Wa * Wb
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
                Ub = Q[i, j, idx, 1] * DbInv - QBar.UBar[i]
                Vb = Q[i, j, idx, 2] * DbInv - QBar.VBar[i]
                Wb = Q[i, j, idx, 3] * DbInv - QBar.WBar[i]
                # Calculate velocity correlation tensor components in the z direction
                R13[kk] = R13[kk] + Ua * Ub
                R23[kk] = R23[kk] + Va * Vb
                R33[kk] = R33[kk] + Wa * Wb
            end
        end
    end
    # Divide by number of points to get averages
    R11 = R11 * nPtsInv
    R21 = R21 * nPtsInv
    R31 = R31 * nPtsInv
    R12 = R12 * nPtsInv
    R22 = R22 * nPtsInv
    R32 = R32 * nPtsInv
    R13 = R13 * nPtsInv
    R23 = R23 * nPtsInv
    R33 = R33 * nPtsInv
    # Normalise by Rab(0)
    i0 = Int(ceil(length(rx) / 2))
    R11 = R11 / R11[i0]
    R21 = R21 / R21[i0]
    R31 = R31 / R31[i0]
    j0 = Int(ceil(length(ry) / 2))
    R12 = R12 / R12[j0]
    R22 = R22 / R22[j0]
    R32 = R32 / R32[j0]
    k0 = Int(ceil(length(rz) / 2))
    R13 = R13 / R13[k0]
    R23 = R23 / R23[k0]
    R33 = R33 / R33[k0]

    return R11, R21, R31, R12, R22, R32, R13, R23, R33

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
    i2 = i0 - 1 + findfirst(x -> x < 0.0, R11[i0:end]) - 1
    i1 = findlast(x -> x < 0.0, R11[1:i0]) + 1
    j2 = j0 - 1 + findfirst(x -> x < 0.0, R22[j0:end]) - 1
    j1 = findlast(x -> x < 0.0, R22[1:j0]) + 1
    k2 = k0 - 1 + findfirst(x -> x < 0.0, R33[k0:end]) - 1
    k1 = findlast(x -> x < 0.0, R33[1:k0]) + 1


    Λx = 0.0
    Λyz = 0.0

    return Λx, Λyz

end