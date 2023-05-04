# PostProc: Post processing routines for CFD simulations on rectilinear grids

# Load modules
import PyPlot as plt
import FortranFiles as FFile
import Dates
import WriteVTK
# Load functions
include("src/file_io.jl")

# Read parameter file
grid, input = readSettings("post.par")
# Get first time step
t = input.startTime
timeStep = rpad(string(round(t, digits=8)), 10, "0")
# Load grid
x, y, z = readPlot3DGrid(timeStep, grid.Nx, grid.Ny, grid.Nz)
# Load solution
Q = readPlot3DSolution(timeStep, grid.Nx, grid.Ny, grid.Nz, input.nVars)

# Write solution to a vtr file 
iNVars = input.nVars
# Check that there are 8, 10 or 12 variables
if (iNVars == 8) || (iNVars == 10) || (iNVars == 12)
    # Write solution to a vtr file
    filename = "data/solution_$(rpad(string(round(t, digits=5)), 7, "0")).vtr"
    tStart = report("Writing the solution at time $(rpad(string(round(t, digits=5)), 7, "0"))", 1)
    WriteVTK.vtk_grid(filename, x[:,1,1], y[1,:,1], z[1,1,:]) do vtk
        vtk["MomentumX"] = Q[:, :, :, 1]
        vtk["MomentumY"] = Q[:, :, :, 2]
        vtk["MomentumZ"] = Q[:, :, :, 3]
        vtk["EnergyDensity"] = Q[:, :, :, 4]
        if (iNVars >= 10)
            for v = 1:2
                vtk["DensityMassFraction$(v)"] = Q[:, :, :, 4+v]
            end
        end
        if (iNVars == 12)
            for v = 1:2
                vtk["VolumeFraction$(v)"] = Q[:, :, :, 4+2+v]
            end
        end
        vtk["Density"] = Q[:, :, :, iNVars-3]
        vtk["Gamma"] = Q[:, :, :, iNVars-2]
        vtk["Pressure"] = Q[:, :, :, iNVars-1]
        vtk["Temperature"] = Q[:, :, :, iNVars]
    end
else
    error("Number of variables must be 8, 10 or 12")
end
tEnd = report("Finished writing solution...", 1)
report("Elapsed time: $(tEnd - tStart)")