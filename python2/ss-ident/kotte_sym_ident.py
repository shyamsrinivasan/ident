import numpy as np
from generate_noisy_data import generate_noisy_data
from generate_noisy_data import run_noisy_parameter_perturbation
from kotte_model import establish_kotte_flux_identifiability
from kotte_model import arrange_experimental_data

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
parameter_perturbation = [(14, 0), (14, 4), (14, 9),
                          (11, .1), (11, .5), (11, 1), (11, -.1), (11, -.5),
                          (12, .1), (12, .5), (12, 1), (12, -.1), (12, -.5),
                          (13, .1), (13, .5), (13, 1), (13, -.1), (13, -.5)]
perturbation_options = {'ode_parameters':ode_paramater_values, 'cvode_options':cvode_options}
noisy_ss, noisy_dynamic, perturbed_parameter_values = \
    run_noisy_parameter_perturbation(parameter_perturbation, noisy_initial_ss[0], perturbation_options)

noisy_exp_xss = []
noisy_exp_fss = []
for ss_values in noisy_ss:
    noisy_exp_xss.append(ss_values[0])
    noisy_exp_fss.append(ss_values[1])

# experimental data based on order of inputs for lambdify expressions
exp_flux_index = np.array([0, 3, 2, 4])
# get combinations of experimental datasets
experimental_datasets = \
    arrange_experimental_data(noisy_exp_xss, noisy_exp_fss, perturbed_parameter_values, exp_flux_index)

# identifiability for all kotte fluxes
establish_kotte_flux_identifiability(experimental_datasets[100:103])

# clear workspace (removes all module names and objects)
# import sys
# sys.modules[__name__].__dict__.clear()
