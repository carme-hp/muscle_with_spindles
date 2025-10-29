import sys
import json
import read_structured_vtk

rank_no = int(sys.argv[-2])
n_ranks = int(sys.argv[-1])

# Simulation input
configuration_name ="configuration1"
isometric = False


if isometric:
    scenario_name = "muscle8_isometric_" + configuration_name
else: 
    scenario_name = "muscle8_isotonic_" + configuration_name

# -------------------------------------------------------------------
# Spindles
# -------------------------------------------------------------------
dt_muscle_spindles = 1e-2
output_timestep_neurons = 1   

n_muscle_spindles = 3

import os
input_dir = os.path.join(os.environ.get('OPENDIHU_HOME', '../../../../../'), "examples/electrophysiology/input/")
muscle_spindle_cellml_file = input_dir + "hodgkin_huxley_1952.cellml"

muscle_spindle_mappings = {
  ("parameter", 0):           "membrane/i_Stim",   # stimulation
  ("connectorSlot", 0): "membrane/V",        # voltage
  ("connectorSlot", 1): "membrane/i_Stim",   # stimulation
}
muscle_spindle_parameters_initial_values = [0]    # [i_Stim]

# -------------------------------------------------------------------
# FEM mesh generation from .vts file
# -------------------------------------------------------------------

geometry_name = "muscle8"
vtk_filename = "../"+geometry_name+"/structured_" + geometry_name + ".vtk"
points, bs_x, bs_y, bs_z = read_structured_vtk.read_structured_vtk(vtk_filename)
el_x, el_y, el_z = int((bs_x-1)/2), int((bs_y-1)/2), int((bs_z-1)/2)

meshes = { # create 3D mechanics mesh
    "mesh3D": {
        "nElements":            [el_x, el_y, el_z],
        "nodePositions":        points,
        "logKey":               "mesh3D",
        "inputMeshIsGlobal":    True,
        "nRanks":               1,
    }
}

# -------------------------------------------------------------------
# fiber mesh generation from .json file
# -------------------------------------------------------------------

fiber_file = "../"+geometry_name+"/fibers"+geometry_name+".json"
with open(fiber_file,"r") as f:
	fdata = json.load(f)

fiber_idx = 0
for fiber in fdata:
	fdict = fdata[fiber]
	npos = [[fdict[ii]['x'],fdict[ii]['y'],fdict[ii]['z']] for ii in range(len(fdict)) ]
	meshName = "fiber{}".format(fiber_idx)
	meshes[meshName] = {
			"nElements":		    [len(fdict)-1],
			"nodePositions":	    npos,
			"inputMeshIsGlobal":	True,
			"nRanks":				n_ranks
	}
	fiber_idx += 1
     
n_fibers = fiber_idx

# -------------------------------------------------------------------
# Boundary conditions
# -------------------------------------------------------------------

contraction_dirichlet_bc = {}
bs_x, bs_y, bs_z = el_x*2+1,el_y*2+1, el_z*2+1


# set Dirichlet BC
k = 0
for j in range(bs_y):
  for i in range(bs_x):
    contraction_dirichlet_bc[k*bs_x*bs_y + j*bs_x + i] = [None,None,0.0,None,None,None]

if isometric:
    k = bs_z-1
    for j in range(bs_y):
        for i in range(bs_x):
            contraction_dirichlet_bc[k*bs_x*bs_y + j*bs_x + i] = [None,None,0.0,None,None,None]
       
# set Neumann BC
k = el_z-1
contraction_neumann_bc = [{"element": k*el_x*el_y + j*el_x + i, "constantVector": [0, 0, 0], "face": "2+"} for j in range(el_y) for i in range(el_x)]

# -------------------------------------------------------------------
# fiber setup
# -------------------------------------------------------------------

specific_states_call_enable_begin = 1.0  # time of first fiber activation
specific_states_call_frequency = 1e-3    # frequency of fiber activation