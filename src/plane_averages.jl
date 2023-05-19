# Functions for calculating and manipulating plane averages

# Function to calculate y-z plane averages
function getPlaneAverages(x::SubArray{Float32,1}, Q::Array{Float32,4}, Nx::Int64, Ny::Int64, Nz::Int64, nVars::Int64, thermo::thermodynamicProperties)
    # Reporting
    tStart = report("Calculating plane averages", 1)
    # Extract thermodynamic properties
    μ = thermo.μ
    W = thermo.W
    # Calculate specific gas constants for each species
    R = 8314.46261815324 ./ W
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
    @inbounds begin
        @batch for k = 1:Nz
            for j = 1:Ny
                @simd for i = 1:Nx
                    rhoBar[i] += Q[i, j, k, nVars-3]
                    UBar[i] += Q[i, j, k, 1]
                    if (nVars >= 10)
                        Y1Bar[i] += Q[i, j, k, 5]
                        muBar[i] += 1.0 / (Q[i, j, k, 5] / μ[1] + (1.0 - Q[i, j, k, 5]) / μ[2])
                        nuBar[i] += 1.0 / (Q[i, j, k, nVars-3] * (Q[i, j, k, 5] / μ[1] + (1.0 - Q[i, j, k, 5]) / μ[2]))
                    end
                    if (nVars == 12)
                        Z1Bar[i] += Q[i, j, k, 7]
                        Z1Z2Bar[i] += Q[i, j, k, 7] * Q[i, j, k, 8]
                    elseif (nVars == 10)
                        Z1 = volumeFraction(Q[i, j, k, 5], R)
                        Z1Bar[i] += Z1
                        Z1Z2Bar[i] += Z1 * (1.0 - Z1)
                    end
                end
            end
        end
    end
    # Divide by number of points to get averages
    @inbounds begin
        @turbo for i = 1:Nx
            rhoBar[i] *= nPtsInv
            UBar[i] *= nPtsInv
            Y1Bar[i] *= nPtsInv
            Z1Bar[i] *= nPtsInv
            Z1Z2Bar[i] *= nPtsInv
            muBar[i] *= nPtsInv
            nuBar[i] *= nPtsInv
        end
    end
    # Package into a struct
    QBar = planeAverage(x, rhoBar, UBar, VBar, WBar, Y1Bar, Z1Bar, Z1Z2Bar, muBar, nuBar)
    # Reporting 
    tEnd = report("Finished calculating plane averages...", 1)
    report("Elapsed time: $(tEnd - tStart)")
    return QBar
end

# Function to convert to primitive variables
function convertSolution!(Q::Array{Float32,4}, Nx::Int64, Ny::Int64, Nz::Int64, nVars::Int64)
    # Reporting
    tStart = report("Converting to primitive variables", 1)
    # Loop over all cells
    @inbounds begin
        for k = 1:Nz
            for j = 1:Ny
                @simd for i = 1:Nx
                    rhoInv = 1.0 / Q[i, j, k, nVars-3]
                    Q[i, j, k, 1] *= rhoInv
                    Q[i, j, k, 2] *= rhoInv
                    Q[i, j, k, 3] *= rhoInv
                    if (nVars >= 10)
                        Q[i, j, k, 5] *= rhoInv
                        Q[i, j, k, 6] *= rhoInv
                    end
                end
            end
        end
    end
    # Reporting 
    tEnd = report("Finished converting to primitive variables...", 1)
    report("Elapsed time: $(tEnd - tStart)")
end

# Function to calculate volume fraction from mass fraction (note: assumes species 1)
function volumeFraction(Y1::Float32, R::PtrArray{Float64,1})
    # Calculate denominator
    sumRY = R[1] * Y1 + R[2] * (1.0 - Y1)
    # Calculate volume fraction
    return R[1] * Y1 / sumRY
end