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
v1max_no_enzyme, k1ac_no_enzyme, k1cat_enzyme, k1ac_enzyme = flux_1_ident_expression(experimental_data)
print('{0:10}'.format('Parameter'), '{0:10}'.format('Numerator'),
      '{0:10}'.format('Denominator'), '{0:10}'.format('Value'))
print('{0:10}'.format('V1max:'), '{0:.5f}'.format(v1max_no_enzyme[0]),
      '{0:.5f}'.format(v1max_no_enzyme[1]), '{0:.5f}'.format(v1max_no_enzyme[2]))
print('{0:10}'.format('K1ac:'), '{0:.5f}'.format(k1ac_no_enzyme[0]),
      '{0:.5f}'.format(k1ac_no_enzyme[1]), '{0:.5f}'.format(k1ac_no_enzyme[2]))
print('{0:10}'.format('k1cat:'), '{0:.5f}'.format(k1cat_enzyme[0]),
      '{0:.5f}'.format(k1cat_enzyme[1]), '{0:.5f}'.format(k1cat_enzyme[2]))
print('{0:10}'.format('K1ac:'), '{0:.5f}'.format(k1ac_enzyme[0]),
      '{0:.5f}'.format(k1ac_enzyme[1]), '{0:.5f}'.format(k1ac_enzyme[2]))

# identifiability value for v2
v2max, k2pep = flux_2_ident_expression(experimental_data)
print('{0:10}'.format('V2max:'), '{0:.5f}'.rjust(10).format(v2max[0]),
      '{0:.5f}'.rjust(10).format(v2max[1]), '{0:.5f}'.rjust(10).format(v2max[2]))
print('{0:10}'.format('K2pep:'), '{0:.5f}'.rjust(10).format(k2pep[0]),
      '{0:.5f}'.rjust(10).format(k2pep[1]), '{0:.5f}'.rjust(10).format(k2pep[2]))

# identifiability value for v3
v3max_1, k3fdp_1, k3pep_1, v3max_2, k3fdp_2, k3pep_2 = flux_3_ident_expression(experimental_data)
print('{0:10}'.format('V3max 1:'), '{0:.5f}'.format(v3max_1[0]),
      '{0:.5f}'.format(v3max_1[1]), '{0:.5f}'.format(v3max_1[2]))
print('{0:10}'.format('V3max 2:'), '{0:.5f}'.format(v3max_2[0]),
      '{0:.5f}'.format(v3max_2[1]), '{0:.5f}'.format(v3max_2[2]))
print('{0:10}'.format('K3fdp 1:'), '{0:.5f}'.format(k3fdp_1[0]),
      '{0:.5f}'.format(k3fdp_1[1]), '{0:.5f}'.format(k3fdp_1[2]))
print('{0:10}'.format('K3fdp 2:'), '{0:.5f}'.format(k3fdp_2[0]),
      '{0:.5f}'.format(k3fdp_2[1]), '{0:.5f}'.format(k3fdp_2[2]))
print('{0:10}'.format('K3pep 1:'), '{0:.5f}'.format(k3pep_1[0]),
      '{0:.5f}'.format(k3pep_1[1]), '{0:.5f}'.format(k3pep_1[2]))
print('{0:10}'.format('K3pep 2:'), '{0:.5f}'.format(k3pep_2[0]),
      '{0:.5f}'.format(k3pep_2[1]), '{0:.5f}'.format(k3pep_2[2]))
