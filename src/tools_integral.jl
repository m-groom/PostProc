# Tools for use in functions from integral_quantities.jl

# Trapezoidal rule integration (x is a vector of nodes, y is a vector of values at nodes)
function trapz(x::Array{Float64,1}, y::Array{Float64,1})
    return sum(0.5 .* (x[2:end] .- x[1:end-1]) .* (y[2:end] .+ y[1:end-1]))
end

# Midpoint rule integration (x is a vector of nodes, y is a vector of values at cell centres)
function midpoint(x::Array{Float64,1}, y::Array{Float64,1})
    return sum((x[2:end] .- x[1:end-1]) .* y[1:end])
end

# Function to compute a first order derivative in the non-homogeneous direction
function dUdX(x::Array{Float32,1}, U::Array{Float32,1}, i::Int64)
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
function dUdY(y::Array{Float32,1}, U::Array{Float32,1}, j::Int64)
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