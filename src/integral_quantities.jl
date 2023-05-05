# Functions for calculating and manipulating integral quantities

# Top level function to calculate integral quantities
function calcIntegralQuantities(t::Float64, x::Array{Float32,3}, y::Array{Float32,3}, z::Array{Float32,3}, Q::Array{Float32,4}, QBar::planeAverage, grid::rectilinearGrid, nVars::Int64, x0::Float64)
    # Calculate velocity correlation tensor
    R11, R22, R33 = calcVelocityCorrelation(x, y, z, Q, QBar, grid, nVars)

end

# Function to calculate velocity correlation tensor
function calcVelocityCorrelation(x::Array{Float32,3}, y::Array{Float32,3}, z::Array{Float32,3}, Q::Array{Float32,4}, QBar::planeAverage, grid::rectilinearGrid, nVars::Int64)
    # Initialise sums (currently just for r=0)
    R11 = zeros(Float64, grid.Nx) #, grid.Nx)
    R22 = zeros(Float64, grid.Nx) #, grid.Ny)
    R33 = zeros(Float64, grid.Nx) #, grid.Nz)
    nPtsInv = 1.0 / (grid.Ny * grid.Nz)
    # Find index where x = xL
    iL = argmin(abs.(x[:, 1, 1] .- grid.xL))
    # Find index where x = xR
    iR = argmin(abs.(x[:, 1, 1] .- grid.xR))
    # Loop over all cells
    for k = 1:grid.Nz
        for j = 1:grid.Ny
            for i = 1:grid.Nx
                # # Get separation directions
                # rx = x[iL:iR, j, k]
                # ry = y[i, :, k]
                # rz = z[i, j, :]
                # Calculate Rii
                DInv = 1.0 / Q[i, j, k, nVars-3]
                u = Q[i, j, k, 1] * DInv
                v = Q[i, j, k, 2] * DInv
                w = Q[i, j, k, 3] * DInv
                R11[i] = R11[i] + (u - QBar.UBar[i])^2
                R22[i] = R22[i] + (v - QBar.VBar[i])^2
                R33[i] = R33[i] + (w - QBar.WBar[i])^2
            end
        end
    end
    # Divide by number of points to get averages
    R11 = R11 * nPtsInv
    R22 = R22 * nPtsInv
    R33 = R33 * nPtsInv

    return R11, R22, R33

end