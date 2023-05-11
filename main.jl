# PostProc: Post processing routines for CFD simulations on rectilinear grids

# Load modules
import PyPlot as plt
import FortranFiles as FFile
import Dates
import WriteVTK
import FFTW
using Printf
# Load functions
include("src/structs.jl")
include("src/file_io.jl")
include("src/plane_averages.jl")
include("src/integral_quantities.jl")
include("src/tools_integral.jl")
include("src/spectral_quantities.jl")
include("src/tools_spectral.jl")

# Reporting
t1 = report("Starting post-processing", 1)
# Read parameter file
grid, input, thermo, x0, dataDir = readSettings("post.par")
# Get first time step
t = input.startTime
timeStep = rpad(string(round(t, digits=8)), 10, "0")
# Load grid
x, y, z = readPlot3DGrid(timeStep, grid.Nx, grid.Ny, grid.Nz, dataDir)

# Loop over all time steps
for n = 1:input.nFiles
    # Get time step
    global t = input.startTime + (n - 1) * input.Δt
    global timeStep = rpad(string(round(t, digits=8)), 10, "0")
    # Load solution
    Q = readPlot3DSolution(timeStep, grid.Nx, grid.Ny, grid.Nz, input.nVars, dataDir)

    # Write out full solution
    writeSolution(t, x, y, z, Q, input.nVars, dataDir)

    # Write out slice
    writeSlice(t, x, y, z, Q, input.nVars, "xy", x0, dataDir)

    # Calculate plane averages
    QBar = getPlaneAverages(x, Q, grid.Nx, grid.Ny, grid.Nz, input.nVars, thermo)

    # Write plane averages
    writePlaneAverages(t, QBar, dataDir)

    # Calculate integral quantities
    calcIntegralQuantities(t, x, y, z, Q, QBar, grid, input.nVars, x0, dataDir)

    # Calculate spectral quantities
    calcSpectralQuantities(t, x, y, z, Q, QBar, grid, input.nVars, x0, dataDir)
end

# Reporting
t2 = report("Finished post-processing...", 1)
report("Total time: $(t2 - t1)")