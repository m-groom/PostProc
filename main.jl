# PostProc: Post processing routines for CFD simulations on rectilinear grids

# Load modules
import PyPlot as plt
import FortranFiles as FFile
import Dates
import WriteVTK
using Printf
# Load functions
include("src/structs.jl")
include("src/file_io.jl")
include("src/plane_averages.jl")
include("src/integral_quantities.jl")
include("src/tools_integral.jl")

# Read parameter file
grid, input, x0, μ = readSettings("post.par")
# Get first time step
t = input.startTime
timeStep = rpad(string(round(t, digits=8)), 10, "0")
# Load grid
x, y, z = readPlot3DGrid(timeStep, grid.Nx, grid.Ny, grid.Nz)

# TODO: wrap this in a loop to read all time steps
t = 0.0
timeStep = rpad(string(round(t, digits=8)), 10, "0")
# Load solution
Q = readPlot3DSolution(timeStep, grid.Nx, grid.Ny, grid.Nz, input.nVars)

# Write out full solution
writeSolution(t, x, y, z, Q, input.nVars)

# Write out slice
writeSlice(t, x, y, z, Q, input.nVars, "xy", x0)

# Calculate plane averages
QBar = getPlaneAverages(x, Q, grid.Nx, grid.Ny, grid.Nz, input.nVars, μ)

# Write plane averages
writePlaneAverages(t, QBar)

# Calculate integral quantities
calcIntegralQuantities(t, x, y, z, Q, QBar, grid, input.nVars, x0)