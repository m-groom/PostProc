# File for storing structs used in PostProc

# Struct to store settings for a rectilinear grid
struct rectilinearGrid
    Nx::Integer
    Ny::Integer
    Nz::Integer
    xL::AbstractFloat
    xR::AbstractFloat
    yL::AbstractFloat
    yR::AbstractFloat
    zL::AbstractFloat
    zR::AbstractFloat
end

# Struct to store input file settings
struct inputSettings
    nVars::Integer
    startTime::AbstractFloat
    Δt::AbstractFloat
    nFiles::Integer
end

# Struct to store thermodynamic properties
struct thermodynamicProperties
    μ::AbstractArray
    W::AbstractArray
end

# Struct to store plane averages
struct planeAverage
    x::AbstractArray
    rhoBar::AbstractArray
    UBar::AbstractArray
    VBar::AbstractArray
    WBar::AbstractArray
    Y1Bar::AbstractArray
    Z1Bar::AbstractArray
    Z1Z2Bar::AbstractArray
    muBar::AbstractArray
    nuBar::AbstractArray
end