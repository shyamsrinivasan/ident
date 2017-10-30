# clear workspace (removes all module names and objects)
import sys
sys.modules[__name__].__dict__.clear()
import numpy as np
from generate_noisy_data import generate_noisy_data
from generate_noisy_data import run_noisy_parameter_perturbation
from kotte_model import flux_1_ident_expression
from kotte_model import flux_2_ident_expression
from kotte_model import flux_3_ident_expression

# generate noisy experimental data for testing identifiability
y0 = np.array([5, 1, 1])
all_options_exp_1 = []
all_options_exp_2 = []
all_options_exp_3 = []
# default parameter values
cvode_options = ('Newton', 'Adams', 1e-6, 1e-6, 100)
ode_paramater_values = np.array([.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, .1])

# get initial noisy system steady state
initial_options = (cvode_options, ode_paramater_values)
noisy_initial_ss, _, _, _ = generate_noisy_data(y0, initial_options, 1)

# all parameter perturbations
parameter_perturbation = [(14, 4), (14, 9)]
perturbation_options = {'ode_parameters':ode_paramater_values, 'cvode_options':cvode_options}
noisy_ss, noisy_dynamic, perturbed_parameter_values = run_noisy_parameter_perturbation(parameter_perturbation, noisy_initial_ss[0], perturbation_options)

noisy_exp_xss = []
noisy_exp_fss = []
for index, ss_values in enumerate(noisy_ss):
    noisy_exp_xss.append(ss_values[0])
    noisy_exp_fss.append(ss_values[1])

# experimental data based on order of inputs for lambdify expressions
exp_flux_index = np.array([0, 3, 2, 4])
experimental_data = np.hstack((ode_paramater_values[-1], noisy_initial_ss[0], noisy_initial_ss[1][exp_flux_index],
                               perturbed_parameter_values[0][-1], noisy_exp_xss[0], noisy_exp_fss[0][exp_flux_index],
                               perturbed_parameter_values[1][-1], noisy_exp_xss[1], noisy_exp_fss[1][exp_flux_index]))

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
