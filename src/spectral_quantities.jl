# Functions for calculating and manipulating spectral quantities

# Top level function to calculate spectral quantities
function calcSpectralQuantities(t::Float64, x::SubArray{Float32,1}, Q::Array{Float32,4}, QBar::planeAverage, grid::rectilinearGrid, dataDir::String)
    # Calculate radial power spectra
    tStart = report("Calculating radial power spectra", 1)
    Eyz, κ = calcPowerSpectra(x, Q, QBar, grid, t, dataDir)
    tEnd = report("Finished calculating radial power spectra...", 1)
    report("Elapsed time: $(tEnd - tStart)")
    # Calculate integral length
    tStart = report("Calculating integral length", 1)
    calcIntegralLength(κ, Eyz, x, grid, t, dataDir)
    tEnd = report("Finished calculating integral length...", 1)
    report("Elapsed time: $(tEnd - tStart)")
end

# Function to calculate radial power spectra at x = x0 for each velocity component
function calcPowerSpectra(x::SubArray{Float32,1}, Q::Array{Float32,4}, QBar::planeAverage, grid::rectilinearGrid, t::Float64, dataDir::String)
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
    # Normalisation constant
    constant = Δy * Δz * min(Δκy, Δκz) / (8.0 * π^2 * grid.Ny * grid.Nz)
    # Calculatate maximum radial wavenumber
    κMax = Int(ceil(sqrt(2.0) * max(grid.Ny / 2, grid.Nz / 2)))
    # Initialise arrays
    uPrime = zeros(Float64, grid.Ny, grid.Nz) # u'
    vPrime = zeros(Float64, grid.Ny, grid.Nz) # v'
    wPrime = zeros(Float64, grid.Ny, grid.Nz) # w'
    E1Dx = zeros(Float64, κMax + 1, grid.Nx) # Energy spectrum for u' at x = x0
    E1Dy = zeros(Float64, κMax + 1, grid.Nx) # Energy spectrum for v' at x = x0
    E1Dz = zeros(Float64, κMax + 1, grid.Nx) # Energy spectrum for w' at x = x0
    κ1D = Float64.(collect(0:κMax)) # radial wavenumber
    # Find index where x = xL
    iL = searchsortedfirst(x, grid.xL)
    # Find index where x = xR
    iR = searchsortedfirst(x, grid.xR)
    # Get fluctuating velocity components
    @inbounds begin
        @batch for i = iL:iR
            @turbo for k = 1:grid.Nz
                for j = 1:grid.Ny
                    uPrime[j, k] = Q[i, j, k, 1] - QBar.UBar[i]
                    vPrime[j, k] = Q[i, j, k, 2] - QBar.VBar[i]
                    wPrime[j, k] = Q[i, j, k, 3] - QBar.WBar[i]
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
            for n in eachindex(κz)
                for m in eachindex(κy)
                    κ = sqrt(κy[m]^2 + κz[n]^2)
                    @simd for p in eachindex(κ1D)
                        κp = p * Δκ
                        if (κp - Δκ / 2 <= κ && κ < κp + Δκ / 2)
                            E1Dx[p+1, i] += E2Dx[m, n] * constant
                            E1Dy[p+1, i] += E2Dy[m, n] * constant
                            E1Dz[p+1, i] += E2Dz[m, n] * constant
                        end
                    end
                end
            end
        end
    end
    # Write energy spectra to file
    writeEnergySpectra(t, κ1D, E1Dx, E1Dy, E1Dz, x, grid, dataDir)
    # Calculate total energy in yz direction
    E1Dyz = @turbo E1Dy .+ E1Dz
    
    return E1Dyz, κ1D

end

# Function to calculate integral length
function calcIntegralLength(κ::Array{Float64,1}, Eyz::Array{Float64,2}, x::SubArray{Float32,1}, grid::rectilinearGrid, t::Float64, dataDir::String)
    # Calculate spacing in κ
    Ly = grid.yR - grid.yL
    Lz = grid.zR - grid.zL
    Δκy = 2.0 * π / Ly
    Δκz = 2.0 * π / Lz
    Δκ = max(Δκy, Δκz)
    # Nyquist wavenumber
    N = Int(max(grid.Ny / 2, grid.Nz / 2)) 
    # Initialise arrays
    Lyz = zeros(Float64, grid.Nx)
    # Find index where x = xL
    iL = searchsortedfirst(x, grid.xL)
    # Find index where x = xR
    iR = searchsortedfirst(x, grid.xR)
    # Store 3π/4 as a constant
    threePiOnFour = 3.0 * π / 4.0
    # Loop over x
    @inbounds begin
        @batch for i = iL:iR
            # Calculate integral length in yz plane
            Lyz[i] = threePiOnFour * integrateEonK(@view(Eyz[:, i]), κ, Δκ, N) / integrate(@view(Eyz[:, i]), Δκ, N)
        end
    end
    # Write integral length to file
    writeIntegralLength(t, Lyz, x, dataDir)

end