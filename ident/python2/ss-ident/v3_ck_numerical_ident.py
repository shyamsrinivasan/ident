 import numpy as np
from generate_expdata import generate_expdata
from simulate_data import arrange_experimental_data_numerical
from numerical_ident import generate_random_initial_conditions
from numerical_ident import solve_multiple_initial_conditions
from plot_ident_results import plot_numerical_parameter_estimates


# generate noisy experimental data for testing identifiability
y0 = np.array([5, 1, 1])
# default parameter values
cvode_options = ('Newton', 'Adams', 1e-10, 1e-10, 200)
ode_parameter_values = {"K1ac": np.array([.1]),
                        "K3fdp": np.array([.1]),
                        "L3fdp": np.array([4e6]),
                        "K3pep": np.array([.1]),
                        "K2pep": np.array([.3]),
                        "vemax": np.array([1.1]),
                        "Kefdp": np.array([.45]),
                        "ne": np.array([2]),
                        "d": np.array([.25]),
                        "V4max": np.array([.2]),
                        "k1cat": np.array([1]),
                        "V3max": np.array([1]),
                        "V2max": np.array([1]),
                        "ac": np.array([.1])}

# get experimental system steady state data without noise using Convenience Kinetics for v3 (kinetics = 2)
exp_xss, exp_fss, exp_ssid, perturbation_details = \
    generate_expdata(y0, cvode_options, ode_parameter_values, noise=0, kinetics=2, dynamic_plot=0,
                     perturbation_plot=0)

# arrange experimental data to form multiple data sets
exp_flux_index = np.array([0, 3, 2, 4, 1, 5])

# get combination of 3 experiments and perform identifiability on all fluxes that require 3 data sets
print('Practical Identifiability Analysis of v3 with 3 parameters: V3max, K3fdp and K3pep \n')
# choose identifiability functions to test
ident_fun_choice = [0]
# get combinations of experimental datasets
experimental_datasets_3_expts, \
    experiment_choice, combination_choice = arrange_experimental_data_numerical(exp_xss, exp_fss, perturbation_details,
                                                                                experiments_per_set=3,
                                                                                flux_id=exp_flux_index,
                                                                                experiment_choice=[0,
                                                                                                   1, 2, 3, 4, 5,
                                                                                                   6, 7, 8, 9, 10,
                                                                                                   16, 17, 18, 19, 20])
# choose one combination of experimental data and solve nlp for that combination

# NLP solver options
optim_options = {"solver": "ipopt",
                 "opts": {"ipopt.tol": 1e-12}}
# optim_options = {"solver": "sqpmethod",\
#                  "opts": {"qpsol": "qpoases"}}
initial_value = [200, 100, 100, 0, 0, 0]
# randomized_initial_values = generate_random_initial_conditions(initial_value, 10, negative=1)
v3_all_x0_parameter_info = solve_multiple_initial_conditions(all_initial_conditions=[initial_value, [10, 1, 1, 0, 0, 0]],
                                                             experimental_data=experimental_datasets_3_expts[0]["values"],
                                                             chosen_fun=0, optim_options=optim_options,
                                                             number_of_parameters=3, flux_id=3, flux_choice=[3])
plot_numerical_parameter_estimates(v3_all_x0_parameter_info)
print("Run complete\n")
