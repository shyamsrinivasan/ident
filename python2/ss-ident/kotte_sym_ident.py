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
parameter_perturbation = [(14, 4), (14, 9),
                          (12, .1), (12, .5), (12, 1), (12, -.1), (12, -.5),
                          (13, .1), (13, .5), (13, 1), (13, -.1), (13, -.5)]
perturbation_options = {'ode_parameters':ode_paramater_values, 'cvode_options':cvode_options}
noisy_ss, noisy_dynamic, perturbed_parameter_values = \
    run_noisy_parameter_perturbation(parameter_perturbation, noisy_initial_ss[0], perturbation_options)

noisy_exp_xss = []
noisy_exp_fss = []
for index, ss_values in enumerate(noisy_ss):
    noisy_exp_xss.append(ss_values[0])
    noisy_exp_fss.append(ss_values[1])

# experimental data based on order of inputs for lambdify expressions
exp_flux_index = np.array([0, 3, 2, 4])
experimental_data_set_1 = \
    np.hstack((ode_paramater_values[-1], noisy_initial_ss[0], noisy_initial_ss[1][exp_flux_index],
               perturbed_parameter_values[0][-1], noisy_exp_xss[0], noisy_exp_fss[0][exp_flux_index],
               perturbed_parameter_values[1][-1], noisy_exp_xss[1], noisy_exp_fss[1][exp_flux_index]))
experimental_data_set_2 = \
    np.hstack((ode_paramater_values[-1], noisy_initial_ss[0], noisy_initial_ss[1][exp_flux_index],
               perturbed_parameter_values[0][-1], noisy_exp_xss[0], noisy_exp_fss[0][exp_flux_index],
               perturbed_parameter_values[2][-1], noisy_exp_xss[2], noisy_exp_fss[2][exp_flux_index]))
experimental_data_set_3 = \
    np.hstack((ode_paramater_values[-1], noisy_initial_ss[0], noisy_initial_ss[1][exp_flux_index],
               perturbed_parameter_values[0][-1], noisy_exp_xss[0], noisy_exp_fss[0][exp_flux_index],
               perturbed_parameter_values[3][-1], noisy_exp_xss[3], noisy_exp_fss[3][exp_flux_index]))
experimental_data_set_4 = \
    np.hstack((ode_paramater_values[-1], noisy_initial_ss[0], noisy_initial_ss[1][exp_flux_index],
               perturbed_parameter_values[0][-1], noisy_exp_xss[0], noisy_exp_fss[0][exp_flux_index],
               perturbed_parameter_values[4][-1], noisy_exp_xss[4], noisy_exp_fss[4][exp_flux_index]))
experimental_data_set_5 = \
    np.hstack((ode_paramater_values[-1], noisy_initial_ss[0], noisy_initial_ss[1][exp_flux_index],
               perturbed_parameter_values[0][-1], noisy_exp_xss[0], noisy_exp_fss[0][exp_flux_index],
               perturbed_parameter_values[5][-1], noisy_exp_xss[5], noisy_exp_fss[5][exp_flux_index]))
experimental_data_set_6 = \
    np.hstack((ode_paramater_values[-1], noisy_initial_ss[0], noisy_initial_ss[1][exp_flux_index],
               perturbed_parameter_values[0][-1], noisy_exp_xss[0], noisy_exp_fss[0][exp_flux_index],
               perturbed_parameter_values[6][-1], noisy_exp_xss[6], noisy_exp_fss[6][exp_flux_index]))
experimental_data_set_7 = \
    np.hstack((ode_paramater_values[-1], noisy_initial_ss[0], noisy_initial_ss[1][exp_flux_index],
               perturbed_parameter_values[4][-1], noisy_exp_xss[4], noisy_exp_fss[4][exp_flux_index],
               perturbed_parameter_values[6][-1], noisy_exp_xss[6], noisy_exp_fss[6][exp_flux_index]))
experimental_data_set_8 = \
    np.hstack((ode_paramater_values[-1], noisy_initial_ss[0], noisy_initial_ss[1][exp_flux_index],
               perturbed_parameter_values[4][-1], noisy_exp_xss[4], noisy_exp_fss[4][exp_flux_index],
               perturbed_parameter_values[5][-1], noisy_exp_xss[5], noisy_exp_fss[5][exp_flux_index]))
experimental_data_set_9 = \
    np.hstack((ode_paramater_values[-1], noisy_initial_ss[0], noisy_initial_ss[1][exp_flux_index],
               perturbed_parameter_values[5][-1], noisy_exp_xss[5], noisy_exp_fss[5][exp_flux_index],
               perturbed_parameter_values[6][-1], noisy_exp_xss[6], noisy_exp_fss[6][exp_flux_index]))
experimental_data_set_10 = \
    np.hstack((ode_paramater_values[-1], noisy_initial_ss[0], noisy_initial_ss[1][exp_flux_index],
               perturbed_parameter_values[7][-1], noisy_exp_xss[7], noisy_exp_fss[7][exp_flux_index],
               perturbed_parameter_values[8][-1], noisy_exp_xss[8], noisy_exp_fss[8][exp_flux_index]))
experimental_data_set_11 = \
    np.hstack((ode_paramater_values[-1], noisy_initial_ss[0], noisy_initial_ss[1][exp_flux_index],
               perturbed_parameter_values[7][-1], noisy_exp_xss[7], noisy_exp_fss[7][exp_flux_index],
               perturbed_parameter_values[9][-1], noisy_exp_xss[9], noisy_exp_fss[9][exp_flux_index]))
experimental_data_set_12 = \
    np.hstack((ode_paramater_values[-1], noisy_initial_ss[0], noisy_initial_ss[1][exp_flux_index],
               perturbed_parameter_values[9][-1], noisy_exp_xss[9], noisy_exp_fss[9][exp_flux_index],
               perturbed_parameter_values[10][-1], noisy_exp_xss[10], noisy_exp_fss[10][exp_flux_index]))
experimental_data_set_13 = \
    np.hstack((ode_paramater_values[-1], noisy_initial_ss[0], noisy_initial_ss[1][exp_flux_index],
               perturbed_parameter_values[5][-1], noisy_exp_xss[5], noisy_exp_fss[5][exp_flux_index],
               perturbed_parameter_values[10][-1], noisy_exp_xss[10], noisy_exp_fss[10][exp_flux_index]))
experimental_data_set_14 = \
    np.hstack((ode_paramater_values[-1], noisy_initial_ss[0], noisy_initial_ss[1][exp_flux_index],
               perturbed_parameter_values[4][-1], noisy_exp_xss[4], noisy_exp_fss[4][exp_flux_index],
               perturbed_parameter_values[10][-1], noisy_exp_xss[10], noisy_exp_fss[10][exp_flux_index]))
experimental_data_set_15 = \
    np.hstack((ode_paramater_values[-1], noisy_initial_ss[0], noisy_initial_ss[1][exp_flux_index],
               perturbed_parameter_values[3][-1], noisy_exp_xss[3], noisy_exp_fss[3][exp_flux_index],
               perturbed_parameter_values[8][-1], noisy_exp_xss[8], noisy_exp_fss[8][exp_flux_index]))

# identifiability for all kotte fluxes
establish_kotte_flux_identifiability([experimental_data_set_1, experimental_data_set_2, experimental_data_set_3,
                                      experimental_data_set_4, experimental_data_set_5, experimental_data_set_6,
                                      experimental_data_set_7, experimental_data_set_8, experimental_data_set_9,
                                      experimental_data_set_10, experimental_data_set_11, experimental_data_set_12,
                                      experimental_data_set_13, experimental_data_set_14, experimental_data_set_15])

# clear workspace (removes all module names and objects)
# import sys
# sys.modules[__name__].__dict__.clear()
