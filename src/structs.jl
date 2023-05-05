# File for storing structs used in PostProc

# Struct to store settings for a rectilinear grid
struct rectilinearGrid
    Nx::Int64
    Ny::Int64
    Nz::Int64
    xL::Float64
    xR::Float64
    yL::Float64
    yR::Float64
    zL::Float64
    zR::Float64
end

# Struct to store input file settings
struct inputSettings
    nVars::Int
    startTime::Float64
    ΔT::Float64
    nFiles::Int
end

# Struct to store plane averages
struct planeAverage
    x::Array{Float64,1}
    DBar::Array{Float64,1}
    UBar::Array{Float64,1}
    VBar::Array{Float64,1}
    WBar::Array{Float64,1}
    Y1Bar::Array{Float64,1}
    Z1Bar::Array{Float64,1}
    Z1Z2Bar::Array{Float64,1}
end