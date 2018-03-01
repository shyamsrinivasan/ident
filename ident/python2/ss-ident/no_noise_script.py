import numpy as np
from generate_expdata import generate_expdata
from simulate_data import arrange_experimental_data
from kotte_model import flux_ident_2_data_combination
from kotte_model import flux_ident_3_data_combination
from process_ident_data import process_info_sample
from plot_ident_results import data_utility_plot
from plot_ident_results import plot_parameter_values
from plot_ident_results import parameter_identifibaility_plot
from plot_ident_results import parameter_experiment_info_plot
from plot_ident_results import parameter_experiment_info_spider


# generate noisy experimental data for testing identifiability
y0 = np.array([5, 1, 1])
# default parameter values
cvode_options = ('Newton', 'Adams', 1e-10, 1e-10, 200)
# ode_parameter_values = np.array([.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, .1])
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

# get experimental system steady state data without noise
exp_xss, exp_fss, exp_ssid, perturbation_details = \
    generate_expdata(y0, cvode_options, ode_parameter_values, noise=0, kinetics=2, dynamic_plot=0,
                     perturbation_plot=0)

# arrange experimental data to form multiple data sets
exp_flux_index = np.array([0, 3, 2, 4, 1, 5])

# get combination of 2 experiments and perform identifiability on all fluxes that require 2 data sets
print('Practical Identifiability Analysis of fluxes with 2 parameters \n')
# choose which identifiability functions to test
ident_fun_choice_2 = [0, 1]
# get combinations of experimental datasets
experimental_datasets_2_expts, \
    experiment_choice_2, combination_choice_2 = arrange_experimental_data(exp_xss, exp_fss, perturbation_details,
                                                                          experiments_per_set=2, flux_id=exp_flux_index,
                                                                          experiment_choice=[0, 1, 2, 8, 9, 10, 11, 12])
ident_details_2 = flux_ident_2_data_combination(experimental_datasets_2_expts, choose=combination_choice_2,
                                                flux_ids=[1, 2], flux_choice=[2, 0], ident_fun_choice=ident_fun_choice_2)
print('Identifiability analysis for fluxes with 2 parameters complete.\n')
# data processing
experiment_type_indices = [[0], [1, 2], [3, 4, 5, 6, 7], [8, 9, 10, 11, 12], [13, 14, 15, 16, 17]]
data_list_2, max_parameter_2, true_value_2, experiment_info_2, \
    combined_data_list_2, combined_max_parameter_2, combined_true_value_2, \
    combined_experiment_info_2 = process_info_sample(ident_details_2,
                                                     experimental_datasets_2_expts,
                                                     experiment_type_indices,
                                                     combine_fluxes=0,
                                                     ident_fun_choice=ident_fun_choice_2)

# plot parameter identifibaility for all fluxes using 2 data combinations
parameter_identifibaility_plot(max_parameter_2)
# plot experiment type in each position based on all parameter
# identifiable data combinations for each parameter
# parameter_experiment_info_plot(experiment_info_2)
parameter_experiment_info_spider(experiment_info_2)
# plot true parameter values and determined parameter values
plot_parameter_values(true_value_2)

# plot utility of data sets (number of data sets identifying n, n-1, n-2, ...., 1, 0 parameters
data_utility_plot(data_list_2)

# plot combined data utility
# data_utility_plot(combined_data_list_2)

# get combination of 3 experiments and perform identifiability on all fluxes that require 3 data sets
print('Practical Identifiability Analysis of fluxes with 3 parameters \n')
# choose identifiability functions to test
ident_fun_choice_3 = [0]
# get combinations of experimental datasets
experimental_datasets_3_expts, \
    experiment_choice_3, combination_choice_3 = arrange_experimental_data(exp_xss, exp_fss, perturbation_details,
                                                                          experiments_per_set=3, flux_id=exp_flux_index,
                                                                          experiment_choice=[0, 1, 2, 3, 4, 5, 6, 7,
                                                                                             13, 14, 15, 16, 17])
ident_details_3 = flux_ident_3_data_combination(experimental_datasets_3_expts, choose=combination_choice_3,
                                                flux_ids=[3], flux_choice=[1], ident_fun_choice=ident_fun_choice_3)
print('Identifiability analysis for fluxes with 2 parameters complete.\n')
# data processing
data_list_3, max_parameter_3, true_value_3, experiment_info_3, \
    combined_data_list_3, combined_max_parameter_3, combined_true_value_3, \
    combined_experiment_info_3 = process_info_sample(ident_details_3,
                                                     experimental_datasets_3_expts,
                                                     experiment_type_indices,
                                                     ident_fun_choice=ident_fun_choice_3)

# plot
# plot parameter identifibaility for all fluxes using 3 data combinations
parameter_identifibaility_plot(max_parameter_3)
# plot experiment type in each position based on all parameter
# identifiable data combinations for each parameter
# parameter_experiment_info_plot(experiment_info_3)
parameter_experiment_info_spider(experiment_info_3)
# plot true parameter values and determined parameter values
plot_parameter_values(true_value_3)

# plot utility of data sets (number of data sets identifying n, n-1, n-2, ...., 1, 0 parameters
data_utility_plot(data_list_3)


# identifiability for all kotte fluxes
# ident_details = establish_kotte_flux_identifiability(experimental_datasets, choose=choose)

# plot results
# plot different experiment types identifying each parameter
# parameter_choice = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

print("Run complete\n")
