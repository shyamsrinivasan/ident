from sympy import *

# generate noisy experimental data for testing identifiability


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