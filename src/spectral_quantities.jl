# Functions for calculating and manipulating spectral quantities

# Top level function to calculate spectral quantities
function calcSpectralQuantities(t::Float64, x::Array{Float32,3}, y::Array{Float32,3}, z::Array{Float32,3}, Q::Array{Float32,4}, QBar::planeAverage, grid::rectilinearGrid, nVars::Int64, x0::Float64)
    # Calcualte radial power spectra
    tStart = report("Calculating radial power spectra", 1)
    Ex, Ey, Ez, κ = calcPowerSpectra(t, x[:, 1, 1], y[1, :, 1], z[1, 1, :], Q, QBar, grid, nVars, x0)
    tEnd = report("Finished calculating radial power spectra...", 1)
    report("Elapsed time: $(tEnd - tStart)")


end

# Function to calculate radial power spectra at x = x0 for each velocity component
function calcPowerSpectra(t::Float64, x::Array{Float32,1}, y::Array{Float32,1}, z::Array{Float32,1}, Q::Array{Float32,4}, QBar::planeAverage, grid::rectilinearGrid, nVars::Int64, x0::Float64)
    # Get wavenumbers
    κy = collect(FFTW.fftshift(FFTW.fftfreq(grid.Ny, grid.Ny)))
    κz = collect(FFTW.fftshift(FFTW.fftfreq(grid.Nz, grid.Nz)))
    # Calculatate maximum radial wavenumber
    κMax = Int(ceil(sqrt(κy[1]^2 + κz[1]^2)))
    # Initialise arrays
    uPrime = zeros(Float64, grid.Ny, grid.Nz) # u'
    vPrime = zeros(Float64, grid.Ny, grid.Nz) # v'
    wPrime = zeros(Float64, grid.Ny, grid.Nz) # w'
    κ1D = collect(0:κMax) # radial wavenumber
    E1Dx = zeros(Float64, κMax + 1) # Energy spectrum for u' at x = x0
    E1Dy = zeros(Float64, κMax + 1) # Energy spectrum for v' at x = x0
    E1Dz = zeros(Float64, κMax + 1) # Energy spectrum for w' at x = x0
    # Get location of interface
    i = argmin(abs.(x .- x0))
    # Loop of cells
    for k = 1:grid.Nz
        for j = 1:grid.Ny
            # Get fluctuating velocity components
            rhoInv = 1.0 / Q[i, j, k, nVars-3]
            uPrime[j, k] = Q[i, j, k, 1] * rhoInv - QBar.UBar[i]
            vPrime[j, k] = Q[i, j, k, 2] * rhoInv - QBar.VBar[i]
            wPrime[j, k] = Q[i, j, k, 3] * rhoInv - QBar.WBar[i]
        end
    end
    # 2D FFT
    uHat = FFTW.fftshift(FFTW.fft(uPrime))
    vHat = FFTW.fftshift(FFTW.fft(vPrime))
    wHat = FFTW.fftshift(FFTW.fft(wPrime))
    E2Dx = conj(uHat) .* uHat
    E2Dy = conj(vHat) .* vHat
    E2Dz = conj(wHat) .* wHat
    # Calculate power spectra
    num = zeros(Float64, κMax + 1)
    for n in eachindex(κz)
        for m in eachindex(κy)
            κ = sqrt(κy[m]^2 + κz[n]^2)
            idx = Int(ceil(κ)) + 1
            E1Dx[idx] += E2Dx[m, n]
            E1Dy[idx] += E2Dy[m, n]
            E1Dz[idx] += E2Dz[m, n]
            num[idx] += 1
        end
    end
    # Normalise
    for κ = 0:κMax
        idx = κ + 1
        if (num[idx] == 0)
            num[idx] = 1
        end
        area = 0.25 * pi * ((κ + 0.5)^2 - (κ - 0.5)^2)
        E1Dx[idx] = E1Dx[idx] * area / num[idx]
        E1Dy[idx] = E1Dy[idx] * area / num[idx]
        E1Dz[idx] = E1Dz[idx] * area / num[idx]
    end
    # TODO print out total energy

    return E1Dx, E1Dy, E1Dz, κ1D

end