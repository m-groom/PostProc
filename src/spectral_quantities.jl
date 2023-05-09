# Functions for calculating and manipulating spectral quantities

# Top level function to calculate spectral quantities
function calcSpectrallQuantities(t::Float64, x::Array{Float32,3}, y::Array{Float32,3}, z::Array{Float32,3}, Q::Array{Float32,4}, QBar::planeAverage, grid::rectilinearGrid, nVars::Int64, x0::Float64)
    # Calcualte radial power spectra
    Ex, Ey, Ez = calcPowerSpectra(t, x, y, z, Q, QBar, grid, nVars, x0)

end

# Function to calculate radial power spectra at x = x0 for each velocity component
function calcPowerSpectra(t::Float64, x::Array{Float32,1}, y::Array{Float32,1}, z::Array{Float32,1}, Q::Array{Float32,4}, QBar::planeAverage, grid::rectilinearGrid, nVars::Int64, x0::Float64)

end