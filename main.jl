# PostProc: Post processing routines for CFD simulations on rectilinear grids

# Load modules
import PyPlot as plt
# Load functions
include("src/file_io.jl")

# Read paraeter file
grid, input = readSettings("post.par")
# # Load data
# x, y, z = loadGrid(filename)
# Q = loadSolution(filename)