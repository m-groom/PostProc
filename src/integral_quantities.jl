# Functions for calculating and manipulating integral quantities

# Top level function to calculate integral quantities
function calcIntegralQuantities(t::AbstractFloat, x::AbstractArray, y::AbstractArray, z::AbstractArray, Q::AbstractArray, QBar::planeAverage, grid::rectilinearGrid, dataDir::AbstractString)
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
    calcLengthScales(t, x, y, z, Q, QBar, grid, dataDir)
    tEnd = report("Finished calculating directional length scales...", 1)
    report("Elapsed time: $(tEnd - tStart)")
end

# Function to calculate velocity correlation tensor at x = x0
function calcVelocityCorrelation(x::AbstractArray, y::AbstractArray, z::AbstractArray, Q::AbstractArray, QBar::planeAverage, grid::rectilinearGrid, t::AbstractFloat, dataDir::AbstractString)
    # Find index where x = xL
    iL = searchsortedfirst(x, grid.xL)
    # Find index where x = xR
    iR = searchsortedfirst(x, grid.xR)
    # Get separation directions
    ry = y .- 0.5 .* y[end]
    rz = z .- 0.5 .* z[end]
    Nx = 2 * (length(@view(x[iL:iR])) - 1)
    rx = collect(range(-grid.xR, stop=grid.xR, length=Nx + 1))
    # Location of r=0
    i0 = Int(ceil(length(rx) / 2))
    j0 = Int(ceil(length(ry) / 2))
    k0 = Int(ceil(length(rz) / 2))
    # Initialise arrays
    R11 = zeros(Float64, length(rx), grid.Nx) # <u'(x)u'(x+rx)>
    R22 = zeros(Float64, length(ry), grid.Nx) # <v'(x)v'(x+ry)>
    R33 = zeros(Float64, length(rz), grid.Nx) # <w'(x)w'(x+rz)>
    # Get number of cells
    Ny = grid.Ny
    Nz = grid.Nz
    nPtsInv = 1.0 / (Ny * Nz)
    # Loop over all cells
    @inbounds begin
        @batch for i = iL:iR
            for k = 1:Nz
                for j = 1:Ny
                    # Calculate fluctuating velocities at (i,j,k)
                    Ua = Q[i, j, k, 1] - QBar.UBar[i]
                    Va = Q[i, j, k, 2] - QBar.VBar[i]
                    Wa = Q[i, j, k, 3] - QBar.WBar[i]
                    # Calculate correlations in the x direction
                    for ii in eachindex(rx)
                        # Get index of point to calculate flucuation at (note: only correct within region of uniform grid spacing)
                        idx = Int(i - 1 + ii - Nx / 2)
                        # Restrict x to be within bounds
                        if (idx < 1)
                            idx = 1
                        elseif (idx > grid.Nx)
                            idx = grid.Nx
                        end
                        # Interpolate velocity at idx
                        Ub = Q[idx, j, k, 1] - QBar.UBar[idx]
                        # Calculate velocity correlation tensor components in the x direction
                        R11[ii, i] += Ua * Ub
                    end
                    # Calculate correlations in the y direction
                    for jj in eachindex(ry)
                        # Get index of point to calculate flucuation at
                        idx = Int(j - 1 + jj - Ny / 2)
                        # Restrict y to be within 0 and 2π
                        if (idx <= 0)
                            idx += Ny
                        elseif (idx > Ny)
                            idx -= Ny
                        end
                        # Interpolate velocity at idx
                        Vb = Q[i, idx, k, 2] - QBar.VBar[i]
                        # Calculate velocity correlation tensor components in the y direction
                        R22[jj, i] += Va * Vb
                    end
                    # Calculate correlations in the z direction
                    for kk in eachindex(rz)
                        # Get index of point to calculate flucuation at
                        idx = Int(k - 1 + kk - Nz / 2)
                        # Restrict z to be within 0 and 2π
                        if (idx <= 0)
                            idx += Nz
                        elseif (idx > Nz)
                            idx -= Nz
                        end
                        # Interpolate velocity at idx
                        Wb = Q[i, j, idx, 3] - QBar.WBar[i]
                        # Calculate velocity correlation tensor components in the z direction
                        R33[kk, i] += Wa * Wb
                    end
                end
            end
            # Divide by number of points to get averages
            R11[:, i] *= nPtsInv
            R22[:, i] *= nPtsInv
            R33[:, i] *= nPtsInv
            # Normalise by Rab(0)
            R11[:, i] /= R11[i0, i]
            R22[:, i] /= R22[j0, i]
            R33[:, i] /= R33[k0, i]
        end
    end
    # Write velocity correlations to file
    writeVelocityCorrelation(t, x, R11, R22, R33, grid, dataDir)

    return R11, R22, R33

end

# Function to calculate (longitudinal) correlation lengths
function calcCorrelationLengths(R11::AbstractArray, R22::AbstractArray, R33::AbstractArray, x::AbstractArray, y::AbstractArray, z::AbstractArray, grid::rectilinearGrid)
    # Find index where x = xL
    iL = searchsortedfirst(x, grid.xL)
    # Find index where x = xR
    iR = searchsortedfirst(x, grid.xR)
    # Get separation directions
    ry = y .- 0.5 .* y[end]
    rz = z .- 0.5 .* z[end]
    Nx = 2 * (length(@view(x[iL:iR])) - 1)
    rx = collect(range(-grid.xR, stop=grid.xR, length=Nx + 1))
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
        for i = iL:iR
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
function calcLengthScales(t::AbstractFloat, x::AbstractArray, y::AbstractArray, z::AbstractArray, Q::AbstractArray, QBar::planeAverage, grid::rectilinearGrid, dataDir::AbstractString)
    # Initialise arrays
    R11 = zeros(Float64, grid.Nx) # <u'u'>
    R22 = zeros(Float64, grid.Nx) # <v'v'>
    R33 = zeros(Float64, grid.Nx) # <w'w'>
    dudxSquared = zeros(Float64, grid.Nx) # <(du/dx)^2>
    dvdySquared = zeros(Float64, grid.Nx) # <(dv/dy)^2>
    dwdzSquared = zeros(Float64, grid.Nx) # <(dw/dz)^2>
    omegaSquared = zeros(Float64, grid.Nx, 3) # <(omega_i)^2>
    divUSquared = zeros(Float64, grid.Nx) # <(du/dx+dv/dy+dw/dz)^2>
    # Get grid spacings and number of points
    Nx = grid.Nx
    Ny = grid.Ny
    Nz = grid.Nz
    ΔyInv = 0.5 / (y[2] - y[1]) # 1/(2Δy)
    ΔzInv = 0.5 / (z[2] - z[1]) # 1/(2Δz)
    nPtsInv = 1.0 / (Ny * Nz)
    # Loop over all cells
    @inbounds begin
        @batch for i = 1:Nx
            for k = 1:Nz
                for j = 1:Ny
                    uPrime = Q[i, j, k, 1] - QBar.UBar[i]
                    vPrime = Q[i, j, k, 2] - QBar.VBar[i]
                    wPrime = Q[i, j, k, 3] - QBar.WBar[i]
                    # Calculate Reynolds stresses
                    R11[i] += uPrime * uPrime
                    R22[i] += vPrime * vPrime
                    R33[i] += wPrime * wPrime
                    # Compute velocity derivatives in x direction
                    if (i == 1) # First order forward difference
                        ΔxInv = 1.0 / (x[2] - x[1])
                        Uii = (Q[2, j, k, 1] - Q[1, j, k, 1]) * ΔxInv
                        Uji = (Q[2, j, k, 2] - Q[1, j, k, 2]) * ΔxInv
                        Uki = (Q[2, j, k, 3] - Q[1, j, k, 3]) * ΔxInv
                    elseif (i == Nx) # First order backward difference
                        ΔxInv = 1.0 / (x[Nx] - x[Nx-1])
                        Uii = (Q[Nx, j, k, 1] - Q[Nx-1, j, k, 1]) * ΔxInv
                        Uji = (Q[Nx, j, k, 2] - Q[Nx-1, j, k, 2]) * ΔxInv
                        Uki = (Q[Nx, j, k, 3] - Q[Nx-1, j, k, 3]) * ΔxInv
                    else # Second order central difference
                        ΔxInv = 1.0 / (x[i+1] - x[i-1])
                        Uii = (Q[i+1, j, k, 1] - Q[i-1, j, k, 1]) * ΔxInv # du/dx
                        Uji = (Q[i+1, j, k, 2] - Q[i-1, j, k, 2]) * ΔxInv # dv/dx
                        Uki = (Q[i+1, j, k, 3] - Q[i-1, j, k, 3]) * ΔxInv # dw/dx
                    end
                    # Compute velocity derivatives in y direction
                    if (j == 1) # Periodic
                        Uij = (Q[i, 2, k, 1] - Q[i, Ny, k, 1]) * ΔyInv
                        Ujj = (Q[i, 2, k, 2] - Q[i, Ny, k, 2]) * ΔyInv
                        Ukj = (Q[i, 2, k, 3] - Q[i, Ny, k, 3]) * ΔyInv
                    elseif (j == Ny) # Periodic
                        Uij = (Q[i, 1, k, 1] - Q[i, Ny-1, k, 1]) * ΔyInv
                        Ujj = (Q[i, 1, k, 2] - Q[i, Ny-1, k, 2]) * ΔyInv
                        Ukj = (Q[i, 1, k, 3] - Q[i, Ny-1, k, 3]) * ΔyInv
                    else # Second order central difference
                        Uij = (Q[i, j+1, k, 1] - Q[i, j-1, k, 1]) * ΔyInv # du/dy
                        Ujj = (Q[i, j+1, k, 2] - Q[i, j-1, k, 2]) * ΔyInv # dv/dy
                        Ukj = (Q[i, j+1, k, 3] - Q[i, j-1, k, 3]) * ΔyInv # dw/dy
                    end
                    # Compute velocity derivatives in z direction
                    if (k == 1) # Periodic
                        Uik = (Q[i, j, 2, 1] - Q[i, j, Nz, 1]) * ΔzInv
                        Ujk = (Q[i, j, 2, 2] - Q[i, j, Nz, 2]) * ΔzInv
                        Ukk = (Q[i, j, 2, 3] - Q[i, j, Nz, 3]) * ΔzInv
                    elseif (k == Nz) # Periodic
                        Uik = (Q[i, j, 1, 1] - Q[i, j, Nz-1, 1]) * ΔzInv
                        Ujk = (Q[i, j, 1, 2] - Q[i, j, Nz-1, 2]) * ΔzInv
                        Ukk = (Q[i, j, 1, 3] - Q[i, j, Nz-1, 3]) * ΔzInv
                    else # Second order central difference
                        Uik = (Q[i, j, k+1, 1] - Q[i, j, k-1, 1]) * ΔzInv # du/dz
                        Ujk = (Q[i, j, k+1, 2] - Q[i, j, k-1, 2]) * ΔzInv # dv/dz
                        Ukk = (Q[i, j, k+1, 3] - Q[i, j, k-1, 3]) * ΔzInv # dw/dz
                    end
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
                    divUSquared[i] += divU * divU
                end
            end
            # Divide by number of points to get averages
            R11[i] *= nPtsInv
            R22[i] *= nPtsInv
            R33[i] *= nPtsInv
            dudxSquared[i] *= nPtsInv
            dvdySquared[i] *= nPtsInv
            dwdzSquared[i] *= nPtsInv
            omegaSquared[i, 1] *= nPtsInv
            omegaSquared[i, 2] *= nPtsInv
            omegaSquared[i, 3] *= nPtsInv
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
function calcDissipationRates(QBar::planeAverage, omegaSquared::AbstractArray, divUSquared::AbstractArray, Nx::Integer)
    # Initialise arrays
    εx = zeros(Float64, Nx)
    εy = zeros(Float64, Nx)
    εz = zeros(Float64, Nx)
    # Store 4/9 as a constant
    fourNinths = 4.0 / 9.0
    # Loop over all cells
    @inbounds begin
        @simd for i = 1:Nx
            nuBar = QBar.muBar[i] / QBar.rhoBar[i]
            εx[i] = nuBar * (omegaSquared[i, 1] + fourNinths * divUSquared[i])
            εy[i] = nuBar * (omegaSquared[i, 2] + fourNinths * divUSquared[i])
            εz[i] = nuBar * (omegaSquared[i, 3] + fourNinths * divUSquared[i])
        end
    end
    return εx, εy, εz
end

# Function to calculate Taylor microscales
function calcTaylorMicroscales(R11::AbstractArray, R22::AbstractArray, R33::AbstractArray, dudxSquared::AbstractArray, dvdySquared::AbstractArray, dwdzSquared::AbstractArray, Nx::Integer)
    # Initialise arrays
    λx = zeros(Float64, Nx)
    λy = zeros(Float64, Nx)
    λz = zeros(Float64, Nx)
    # Loop over all cells
    @inbounds begin
        @simd for i = 1:Nx
            λx[i] = sqrt(R11[i] / dudxSquared[i])
            λy[i] = sqrt(R22[i] / dvdySquared[i])
            λz[i] = sqrt(R33[i] / dwdzSquared[i])
        end
    end
    return λx, λy, λz
end

# Function to calculate Kolmogorov microscales
function calcKolmogorovMicroscales(nuBar::AbstractArray, εx::AbstractArray, εy::AbstractArray, εz::AbstractArray, Nx::Integer)
    # Initialise arrays
    ηx = zeros(Float64, Nx)
    ηy = zeros(Float64, Nx)
    ηz = zeros(Float64, Nx)
    # Loop over all cells
    @inbounds begin
        @simd for i = 1:Nx
            nuBarCubed = nuBar[i]^3
            ηx[i] = (nuBarCubed / εx[i])^0.25
            ηy[i] = (nuBarCubed / εy[i])^0.25
            ηz[i] = (nuBarCubed / εz[i])^0.25
        end
    end
    return ηx, ηy, ηz
end