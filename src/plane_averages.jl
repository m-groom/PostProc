# Functions for calculating and manipulating plane averages

# Function to calculate y-z plane averages
function getPlaneAverages(x::Array{Float32,3}, Q::Array{Float32,4}, Nx::Int64, Ny::Int64, Nz::Int64, nVars::Int64)
    # Initialise sums
    DBar = zeros(Float64, Nx)
    UBar = zeros(Float64, Nx)
    VBar = zeros(Float64, Nx)
    WBar = zeros(Float64, Nx)
    Y1Bar = zeros(Float64, Nx)
    Z1Bar = zeros(Float64, Nx)
    Z1Z2Bar = zeros(Float64, Nx)
    nPtsInv = 1.0 / (Ny * Nz)
    # Loop over all cells
    for k = 1:Nz
        for j = 1:Ny
            for i = 1:Nx
                DInv = 1.0 / Q[i, j, k, nVars-3]
                DBar[i] = DBar[i] + Q[i, j, k, nVars-3]
                UBar[i] = UBar[i] + Q[i, j, k, 1] * DInv
                if (nVars >= 10)
                    Y1Bar[i] = Y1Bar[i] + Q[i, j, k, 5] * DInv
                end
                if (nVars == 12)
                    Z1Bar[i] = Z1Bar[i] + Q[i, j, k, 7]
                    Z1Z2Bar[i] = Z1Z2Bar[i] + Q[i, j, k, 7] * Q[i, j, k, 8]
                end
            end
        end
    end
    # Divide by number of points to get averages
    DBar = DBar * nPtsInv
    UBar = UBar * nPtsInv
    Y1Bar = Y1Bar * nPtsInv
    Z1Bar = Z1Bar * nPtsInv
    Z1Z2Bar = Z1Z2Bar * nPtsInv
    # Package into a struct
    QBar = planeAverage(x[:, 1, 1], DBar, UBar, VBar, WBar, Y1Bar, Z1Bar, Z1Z2Bar)
    return QBar
end

# Function to write out plane averages to a space delimited text file
function writePlaneAverages(t::Float64, QBar::planeAverage)
    # Make file name
    filename = "data/planeAverages_$(rpad(string(round(t, digits=5)), 7, "0")).dat"
    tStart = report("Writing plane averages to file $filename", 1)
    # Open file
    f = open(filename, "w")
    # Write header
    write(f, "# x   DBar   UBar   Y1Bar   Z1Bar   Z1Z2Bar\n")
    # Write data in scientific format with 15 digits
    for i = 1:length(QBar.x)-1
        write(f, "$(@sprintf("%.15e", QBar.x[i]))   $(@sprintf("%.15e", QBar.DBar[i]))   $(@sprintf("%.15e", QBar.UBar[i]))   $(@sprintf("%.15e", QBar.Y1Bar[i]))   $(@sprintf("%.15e", QBar.Z1Bar[i]))   $(@sprintf("%.15e", QBar.Z1Z2Bar[i]))\n")
    end
    # Close file
    close(f)
    tEnd = report("Finished writing plane averages", 1)
    report("Elapsed time: $(tEnd - tStart)")
end

# # Function to write out Reynolds stresses to a space delimited text file
# function writeReynoldsStress(t::Float64, R11::Array{Float64,1}, R22::Array{Float64,1}, R33::Array{Float64,1}, x::Array{Float32,1})
#     # Make file name
#     filename = "data/Rii_$(rpad(string(round(t, digits=5)), 7, "0")).dat"
#     # Open file
#     f = open(filename, "w")
#     # Write header
#     write(f, "# x   R11   R22   R33\n")
#     # Write data in scientific format with 15 digits
#     for i = 1:length(x)-1
#         write(f, "$(@sprintf("%.15e", x[i]))   $(@sprintf("%.15e", R11[i]))   $(@sprintf("%.15e", R22[i]))   $(@sprintf("%.15e", R33[i]))\n")
#     end
#     # Close file
#     close(f)
# end