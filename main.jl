# PostProc: Post processing routines for CFD simulations on rectilinear grids

# Load modules
import PyPlot as plt
# Load functions
include("src/file_io.jl")

# Read parameter file
grid, input = readSettings("post.par")
# # Load data
# x, y, z = loadGrid(filename)
# Q = loadSolution(filename)

x = zeros(grid.Nx+1, grid.Ny+1, grid.Nz+1)
y = zeros(grid.Nx+1, grid.Ny+1, grid.Nz+1)
z = zeros(grid.Nx+1, grid.Ny+1, grid.Nz+1)
Q = zeros(grid.Nx+1, grid.Ny+1, grid.Nz+1, input.nVars)

# Call subroutine readPerProcPlot3D(sTimeStep, iNVar, iNX, iNY, iNZ, rX, rY, rZ, rQ) using ccall
ccall((:readperprocplot3d_, "lib/libreadplot3d.so"), Cvoid, (String, Int64, Int64, Int64, Int64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}), "0.00000000", input.nVars, grid.Nx, grid.Ny, grid.Nz, x, y, z, Q)