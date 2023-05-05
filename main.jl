# PostProc: Post processing routines for CFD simulations on rectilinear grids

# Load modules
import PyPlot as plt
import FortranFiles as FFile
import Dates
import WriteVTK
# Load functions
include("src/file_io.jl")
include("src/structs.jl")

# Read parameter file
grid, input, x0 = readSettings("post.par")
# Get first time step
t = input.startTime
timeStep = rpad(string(round(t, digits=8)), 10, "0")
# Load grid
x, y, z = readPlot3DGrid(timeStep, grid.Nx, grid.Ny, grid.Nz)

# TODO: wrap this in a loop to read all time steps

# Load solution
Q = readPlot3DSolution(timeStep, grid.Nx, grid.Ny, grid.Nz, input.nVars)

# Write out full solution
writeSolution(t, x, y, z, Q, input.nVars)

# Write out slice
# writeSlice(t, x, y, z, Q, input.nVars, "xz", x0)