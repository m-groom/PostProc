# Tools for use in functions from integral_quantities.jl

# Trapezoidal rule
function trapz(x::Array{Float64,1}, y::Array{Float64,1})
    return sum(0.5 .* (x[2:end] .- x[1:end-1]) .* (y[2:end] .+ y[1:end-1]))
end

# Function to write correlation lengths to file
function writeCorrelationLengths(t::Float64, Λx::Float64, Λyz::Float64)
    # Get filename
    filename = "data/correlationLengths.dat"
    report("Writing correlation lengths to file $filename")
    # Append to file
    f = open(filename, "a")
    # Write header
    write(f, "# t   Lx   Lyz\n")
    # Write data in scientific format with 15 digits
    write(f, "$(@sprintf("%.15e", t))   $(@sprintf("%.15e", Λx))   $(@sprintf("%.15e", Λyz))\n")
    # Close file
    close(f)
end