# Tools for use in functions from spectral_quantities.jl

# Function to integrate spectra 
function integrate(E::Array{Float64,1}, Δκ::Float64)
    return sum(E)*Δκ
end