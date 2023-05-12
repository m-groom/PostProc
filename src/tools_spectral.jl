# Tools for use in functions from spectral_quantities.jl

# Function to integrate E(κ)
function integrate(E::Array{Float64,1}, Δκ::Float64, N::Int64)
    integral = 0.0
    @inbounds for i in 1:N
        integral += 0.5 * Δκ * (E[i] + E[i+1])
    end
    return integral
end

# Function to integrate E(κ) / κ
function integrateEonK(E::Array{Float64,1}, κ::Array{Float64,1}, Δκ::Float64, N::Int64)
    integral = 0.0
    @inbounds for i in 1:N
        if (κ[i] == 0)
            integral += 0.5 * Δκ * (0.0 + E[i+1] / κ[i+1])
        else
            integral += +0.5 * Δκ * (E[i] / κ[i] + E[i+1] / κ[i+1])
        end
    end
    return integral
end