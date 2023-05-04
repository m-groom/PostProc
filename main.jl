# PostProc: Post processing routines for CFD simulations on rectilinear grids

# Load modules
import PyPlot as plt
import FortranFiles as FFile
import WriteVTK
# Load functions
include("src/file_io.jl")

# Read parameter file
grid, input = readSettings("post.par")
# Load grid
x, y, z = readPlot3DGrid("0.00000000") 
# Load solution
Q = zeros(size(x))                  



# Write x, y, z to a vts file
WriteVTK.vtk_grid("grid.vts", x, y, z) do vtk
    vtk["Q"] = Q
end