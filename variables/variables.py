import sys
import itertools
import numpy as np

if sys.argv[-1] != "BayesOpt.py":
    n_ranks = int(sys.argv[-1])
    rank_no = int(sys.argv[-2])
else:
    n_ranks = 1

scenario_name = "incompressible_mooney_rivlin"

force = 10.0
constant_body_force = None

# Time stepping
dt_3D = 1e-1            # time step of 3D mechanics
dt_splitting = 2e-3     # time step of strang splitting
dt_1D = 2e-3            # time step of 1D fiber diffusion
dt_0D = 1e-3            # time step of 0D cellml problem
end_time = 40.0         # end time of the simulation 
output_interval = dt_3D # time interval between outputs

# Material parameters
pmax = 7.3                                                  # maximum active stress
rho = 10                                                    # density of the muscle
material_parameters = [3.176e-10, 1.813, 1.075e-2, 1.0]     # [c1, c2, b, d]
diffusion_prefactor = 3.828 / (500.0 * 0.58)                # Conductivity / (Am * Cm)

# Meshes
physical_extent = [3.0, 3.0, 12.0]              # extent of muscle
el_x, el_y, el_z = 3, 3, 12                     # number of elements
bs_x, bs_y, bs_z = 2*el_x+1, 2*el_y+1, 2*el_z+1 # quadratic basis functions

fb_x, fb_y = 10, 10         # number of fibers
fb_points = 100             # number of points per fiber
fiber_direction = [0, 0, 1] # direction of fiber in element

def get_fiber_no(fiber_x, fiber_y):
    return fiber_x + fiber_y*fb_x


# Define directory for cellml files
import os
input_dir = os.path.join(os.environ.get('OPENDIHU_HOME', '../../../../../'), "examples/electrophysiology/input/")

# Fiber activation
fiber_distribution_file = input_dir + "MU_fibre_distribution_3780.txt"
firing_times_file = input_dir + "MU_firing_times_always.txt"
specific_states_call_enable_begin = 1.0                     # time of first fiber activation
specific_states_call_frequency = 1e-3                       # frequency of fiber activation
