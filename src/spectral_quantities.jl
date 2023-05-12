# Functions for calculating and manipulating spectral quantities

# Top level function to calculate spectral quantities
function calcSpectralQuantities(t::Float64, x::SubArray{Float32,1}, Q::Array{Float32,4}, QBar::planeAverage, grid::rectilinearGrid, nVars::Int64, x0::Float64, dataDir::String)
    # Calcualte radial power spectra
    tStart = report("Calculating radial power spectra", 1)
    Ex, Ey, Ez, κ, num = calcPowerSpectra(x, Q, QBar, grid, nVars, x0)
    tEnd = report("Finished calculating radial power spectra...", 1)
    report("Elapsed time: $(tEnd - tStart)")
    # Write energy spectra to file
    writeEnergySpectra(t, κ, Ex, Ey, Ez, grid, num, dataDir)
    # Calculate integral length
    Lyz = calcIntegralLength(κ, Ey, Ez, grid)
    # Write integral length to file
    writeIntegralLength(t, Lyz, dataDir)
end

# Function to calculate radial power spectra at x = x0 for each velocity component
function calcPowerSpectra(x::SubArray{Float32,1}, Q::Array{Float32,4}, QBar::planeAverage, grid::rectilinearGrid, nVars::Int64, x0::Float64)
    # Get wavenumbers and spacing
    κy = collect(FFTW.fftshift(FFTW.fftfreq(grid.Ny, grid.Ny)))
    κz = collect(FFTW.fftshift(FFTW.fftfreq(grid.Nz, grid.Nz)))
    Ly = grid.yR - grid.yL
    Lz = grid.zR - grid.zL
    Δκy = 2.0 * π / Ly
    Δκz = 2.0 * π / Lz
    Δκ = max(Δκy, Δκz)
    Δy = Ly / grid.Ny
    Δz = Lz / grid.Nz
    # Calculatate maximum radial wavenumber
    κMax = Int(ceil(sqrt(2.0) * max(grid.Ny / 2, grid.Nz / 2)))
    # Initialise arrays
    uPrime = zeros(Float64, grid.Ny, grid.Nz) # u'
    vPrime = zeros(Float64, grid.Ny, grid.Nz) # v'
    wPrime = zeros(Float64, grid.Ny, grid.Nz) # w'
    E1Dx = zeros(Float64, κMax + 1) # Energy spectrum for u' at x = x0
    E1Dy = zeros(Float64, κMax + 1) # Energy spectrum for v' at x = x0
    E1Dz = zeros(Float64, κMax + 1) # Energy spectrum for w' at x = x0
    κ1D = Float64.(collect(0:κMax)) # radial wavenumber
    # Find index where x = x0
    i = searchsortedfirst(x, x0) - 1
    # Get fluctuating velocity components
    @inbounds for k = 1:grid.Nz
        @inbounds for j = 1:grid.Ny
            rhoInv = 1.0 / Q[i, j, k, nVars-3]
            uPrime[j, k] = Q[i, j, k, 1] * rhoInv - QBar.UBar[i]
            vPrime[j, k] = Q[i, j, k, 2] * rhoInv - QBar.VBar[i]
            wPrime[j, k] = Q[i, j, k, 3] * rhoInv - QBar.WBar[i]
        end
    end
    # 2D FFT to get power spectral density
    uHat = FFTW.fftshift(FFTW.fft(uPrime))
    vHat = FFTW.fftshift(FFTW.fft(vPrime))
    wHat = FFTW.fftshift(FFTW.fft(wPrime))
    E2Dx = real.(uHat .* conj(uHat))
    E2Dy = real.(vHat .* conj(vHat))
    E2Dz = real.(wHat .* conj(wHat))
    # Calculate energy spectra
    constant = Δy * Δz * min(Δκy, Δκz) / (8.0 * π^2 * grid.Ny * grid.Nz)
    num = zeros(Int, κMax + 1)
    @inbounds for n in eachindex(κz)
        @inbounds for m in eachindex(κy)
            κ = sqrt(κy[m]^2 + κz[n]^2)
            @inbounds for p in eachindex(κ1D)
                κp = p * Δκ
                if (κp - Δκ / 2 <= κ && κ < κp + Δκ / 2)
                    E1Dx[p+1] += E2Dx[m, n] * constant
                    E1Dy[p+1] += E2Dy[m, n] * constant
                    E1Dz[p+1] += E2Dz[m, n] * constant
                    num[p+1] += 1
                end
            end
        end
    end

    return E1Dx, E1Dy, E1Dz, κ1D, num

end

# Function to calculate integral length
function calcIntegralLength(κ::Array{Float64,1}, Ey::Array{Float64,1}, Ez::Array{Float64,1}, grid::rectilinearGrid)
    # Calculate spacing in κ
    Ly = grid.yR - grid.yL
    Lz = grid.zR - grid.zL
    Δκy = 2.0 * π / Ly
    Δκz = 2.0 * π / Lz
    Δκ = max(Δκy, Δκz)
    N = Int(max(grid.Ny / 2, grid.Nz / 2)) # Nyquist wavenumber
    # Calculate integral length in yz plane
    Eyz = Ey .+ Ez
    Lyz = (3.0 * π / 4) * integrateEonK(Eyz, κ, Δκ, N) / integrate(Eyz, Δκ, N)

    return Lyz

end