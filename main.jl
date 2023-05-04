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
x, y, z = readPlot3DGrid("0.00000000", grid.Nx, grid.Ny, grid.Nz) 
# Load solution
Q = readPlot3DSolution("0.00000000", grid.Nx, grid.Ny, grid.Nz, input.nVars)        



# Write solution to a vts file
Q0 = zeros(size(x))   
WriteVTK.vtk_grid("grid.vts", x, y, z) do vtk
    vtk["Q"] = Q0
end