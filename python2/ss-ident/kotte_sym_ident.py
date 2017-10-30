import numpy as np
from generate_noisy_data import generate_noisy_data
from generate_noisy_data import run_noisy_parameter_perturbation
from kotte_model import establish_kotte_flux_identifiability

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

# identifiability for all kotte fluxes
establish_kotte_flux_identifiability([experimental_data])

# clear workspace (removes all module names and objects)
# import sys
# sys.modules[__name__].__dict__.clear()
