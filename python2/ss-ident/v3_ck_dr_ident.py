import numpy as np
from sympy import *
from generate_noisy_data import generate_noisy_data

# generate noisy experimental data for testing identifiability
y0 = np.array([5, 1, 1])
all_options_exp_1 = []
all_options_exp_2 = []
all_options_exp_3 = []
cvode_options = ['Newton', 'Adams', 1e-6, 1e-6, 100]

# experiment 1
ode_par_val_experiment_1 = np.array([.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, .1])
all_options_exp_1.append(cvode_options)
all_options_exp_1.append(ode_par_val_experiment_1)
# generate data using MWC Kinetics
_, y_nss_exp1, flux_nss_exp1, _, _, _, _, _, _ = generate_noisy_data(y0, all_options_exp_1, 1)

# experiment 2
ode_par_val_experiment_2 = np.array([.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, .5])
all_options_exp_2.append(cvode_options)
all_options_exp_2.append(ode_par_val_experiment_2)
# generate data using MWC Kinetics
_, y_nss_exp2, flux_nss_exp2, _, _, _, _, _, _ = generate_noisy_data(y0, all_options_exp_2, 1)

# experiment 3
ode_par_val_experiment_3 = np.array([.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, 1])
all_options_exp_3.append(cvode_options)
all_options_exp_3.append(ode_par_val_experiment_3)
# generate data using MWC Kinetics
_, y_nss_exp3, flux_nss_exp3, _, _, _, _, _, _ = generate_noisy_data(y0, all_options_exp_3, 1)

# experimental data based on order of inputs for lambdify expressions
experimental_data = np.hstack((y_nss_exp1[0:2], flux_nss_exp1[2],
                               y_nss_exp2[0:2], flux_nss_exp2[2],
                               y_nss_exp3[0:2], flux_nss_exp3[2]))

# symbolic expression for flux 3
x11, x12, x13, x21, x22, x23, v31, v32, v33 = symbols('x11, x12, x13, x21, x22, x23, v31, v32, v33', positive=True)
variables = [x11, x21, v31, x12, x22, v32, x13, x23, v33]

# use denominator generated from mathematica to test identifiability for all fluxes
# V3max
V3max_sol_1 = -v32*v33*x11*x12*x21 + v32*v33*x11*x13*x21 + v31*v33*x11*x12*x22 - \
              v31*v33*x12*x13*x22 - v31*v32*x11*x13*x23 + v31*v32*x12*x13*x23
v3max_fun_expression = lambdify([variables], V3max_sol_1, "numpy")
print(v3max_fun_expression(experimental_data))

# K3fdp
# K3fdp_sol_1 - dr same as V3max_sol_1

# K3pep
K3pep_sol_1 = 2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 -
                 v31*v32*x13*x21*x23 - v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)
k3pep_fun_expression = lambdify([variables], K3pep_sol_1, "numpy")
print(k3pep_fun_expression(experimental_data))

# K3fdp_sol_2
# K3pep_sol_2
# V3max_sol_2 = -v32*v33*x11*x12*x21 + v32*v33*x11*x13*x21 + v31*v33*x11*x12*x22 - \
#               v31*v33*x12*x13*x22 - v31*v32*x11*x13*x23 + v31*v32*x12*x13*x23
