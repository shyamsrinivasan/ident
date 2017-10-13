from sympy import *
import numpy as np
from generate_noisy_data import generate_noisy_data

# generate noisy experimental data for testing identifiability
y0 = np.array([5, 1, 1])
# generate data using MWC Kinetics
    _, y_noisy_steady_state, flux_noisy_steady_state, _, _, y_steady_state, flux_steady_state = \
        generate_noisy_data(y0, 1)

variables = [x11, x21, v31, x12, x22, v32, x13, x23, v33]

# use denominator generated from mathematica to test identifiability for all fluxes
# K3fdp_sol_1
# K3pep_sol_1
V3max_sol_1 = -v32*v33*x11*x12*x21 + v32*v33*x11*x13*x21 + v31*v33*x11*x12*x22 - \
              v31*v33*x12*x13*x22 - v31*v32*x11*x13*x23 + v31*v32*x12*x13*x23

# K3fdp_sol_2
# K3pep_sol_2
V3max_sol_2 = -v32*v33*x11*x12*x21 + v32*v33*x11*x13*x21 + v31*v33*x11*x12*x22 - \
              v31*v33*x12*x13*x22 - v31*v32*x11*x13*x23 + v31*v32*x12*x13*x23