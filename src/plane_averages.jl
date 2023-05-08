# Functions for calculating and manipulating plane averages

# Function to calculate y-z plane averages
function getPlaneAverages(x::Array{Float32,3}, Q::Array{Float32,4}, Nx::Int64, Ny::Int64, Nz::Int64, nVars::Int64, μ::Array{Float64,1})
    # Initialise sums
    rhoBar = zeros(Float64, Nx)
    UBar = zeros(Float64, Nx)
    VBar = zeros(Float64, Nx)
    WBar = zeros(Float64, Nx)
    Y1Bar = zeros(Float64, Nx)
    Z1Bar = zeros(Float64, Nx)
    Z1Z2Bar = zeros(Float64, Nx)
    muBar = zeros(Float64, Nx)
    nuBar = zeros(Float64, Nx)
    nPtsInv = 1.0 / (Ny * Nz)
    # Loop over all cells
    for k = 1:Nz
        for j = 1:Ny
            for i = 1:Nx
                rhoInv = 1.0 / Q[i, j, k, nVars-3]
                rhoBar[i] = rhoBar[i] + Q[i, j, k, nVars-3]
                UBar[i] = UBar[i] + Q[i, j, k, 1] * rhoInv
                if (nVars >= 10)
                    Y1 = Q[i, j, k, 5] * rhoInv
                    Y1Bar[i] = Y1Bar[i] + Y1
                    muBar[i] = muBar[i] + 1.0 / (Y1 / μ[1] + (1.0 - Y1) / μ[2])
                    nuBar[i] = nuBar[i] + rhoInv / (Y1 / μ[1] + (1.0 - Y1) / μ[2])
                end
                if (nVars == 12)
                    Z1Bar[i] = Z1Bar[i] + Q[i, j, k, 7]
                    Z1Z2Bar[i] = Z1Z2Bar[i] + Q[i, j, k, 7] * Q[i, j, k, 8]
                end
            end
        end
    end
    # Divide by number of points to get averages
    rhoBar = rhoBar * nPtsInv
    UBar = UBar * nPtsInv
    Y1Bar = Y1Bar * nPtsInv
    Z1Bar = Z1Bar * nPtsInv
    Z1Z2Bar = Z1Z2Bar * nPtsInv
    muBar = muBar * nPtsInv
    nuBar = nuBar * nPtsInv
    # Package into a struct
    QBar = planeAverage(x[:, 1, 1], rhoBar, UBar, VBar, WBar, Y1Bar, Z1Bar, Z1Z2Bar, muBar, nuBar)
    return QBar
end