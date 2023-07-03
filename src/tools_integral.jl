# Tools for use in functions from integral_quantities.jl

# Trapezoidal rule integration (x is a vector of nodes, y is a vector of values at nodes)
function trapz(x::AbstractArray, y::AbstractArray)
    integral = 0.0
    @inbounds begin
        @simd for i in firstindex(x):lastindex(x)-1
            integral += 0.5 * (x[i+1] - x[i]) * (y[i+1] + y[i])
        end
    end
    return integral
end

# Midpoint rule integration (x is a vector of nodes, y is a vector of values at cell centres)
function midpoint(x::AbstractArray, y::AbstractArray)
    integral = 0.0
    @inbounds begin
        @simd for i in firstindex(x):lastindex(x)-1
            integral += (x[i+1] - x[i]) * y[i]
        end
    end
    return integral
end

# Function to compute a first order derivative in the non-homogeneous direction
function dUdX(x::AbstractArray, U::AbstractArray, i::Integer)
    # Get number of cells
    N = length(x) - 1
    # Check if we are at the first or last cell
    if i == 1 # First order
        return (U[2] - U[1]) / (x[2] - x[1])
    elseif i == N # First order
        return (U[N] - U[N-1]) / (x[N] - x[N-1])
    else
        return (U[i+1] - U[i-1]) / (x[i+1] - x[i-1])
    end
end

# Function to compute a first order derivative in the homogeneous direction
function dUdY(y::AbstractArray, U::AbstractArray, j::Integer)
    # Get number of cells
    N = length(y) - 1
    # Check if we are at the first or last cell
    if j == 1 # Periodic
        return (U[2] - U[N]) / (y[2] - y[N])
    elseif j == N # Periodic
        return (U[1] - U[N-1]) / (y[1] - y[N-1])
    else
        return (U[j+1] - U[j-1]) / (y[j+1] - y[j-1])
    end
end

# Function to interpolate the x velocity at a given point
function interp1D(x0::AbstractFloat, x::AbstractArray, U::AbstractArray, UBar::AbstractArray)
    # Find the cell containing x0
    i = searchsortedfirst(x, x0) - 1
    # Grid spacings
    ΔxC = x[i+1] - x[i]
    ΔxL = x[i] - x[i-1]
    ΔxR = x[i+2] - x[i+1]
    # Get values at node i
    Um = ((U[i-1] - UBar[i-1]) * ΔxL + (U[i] - UBar[i]) * ΔxC) / (ΔxL + ΔxC)
    # Get values at node i+1
    Up = ((U[i] - UBar[i]) * ΔxC + (U[i+1] - UBar[i+1]) * ΔxR) / (ΔxC + ΔxR)
    # Interpolate
    return Up + (Up - Um) * (x0 - x[i]) / ΔxC
end