import numpy as np
from generate_noisy_data import generate_noisy_data
from kotte_model import flux_1_ident_expression
from kotte_model import flux_2_ident_expression
from kotte_model import flux_3_ident_expression

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
exp_flux_index = np.array([0, 3, 2, 4])
experimental_data = np.hstack((ode_par_val_experiment_1[-1], y_nss_exp1, flux_nss_exp1[exp_flux_index],
                               ode_par_val_experiment_2[-1], y_nss_exp2, flux_nss_exp2[exp_flux_index],
                               ode_par_val_experiment_3[-1], y_nss_exp3, flux_nss_exp3[exp_flux_index]))

# identifiability value for v1
no_enzyme_dr, enzyme_dr = flux_1_ident_expression(experimental_data)
print("V1max Denominator (No enzyme data):", no_enzyme_dr[0])
print("K1ac Denominator (No enzyme data):", no_enzyme_dr[1])
print("k1cat Denominator (w/ enzyme data):", enzyme_dr[0])
print("K1ac Denominator (w/ enzyme data):", enzyme_dr[1])

# identifiability value for v2
v2_no_enzyme_dr = flux_2_ident_expression(experimental_data)
print("V2max Denominator (No enzyme data):", v2_no_enzyme_dr[0])
print("K2pep Denominator (No enzyme data):", v2_no_enzyme_dr[1])

# identifiability value for v3
v3_no_enzyme_dr = flux_3_ident_expression(experimental_data)
print("V3max Denominator (No enzyme data):", v3_no_enzyme_dr[0])
print("K3fdp Denominator (No enzyme data):", v3_no_enzyme_dr[1])
print("K3pep Denominator (No enzyme data):", v3_no_enzyme_dr[2])
