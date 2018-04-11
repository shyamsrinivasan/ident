import numpy as np
from generate_expdata import generate_expdata
from simulate_data import arrange_experimental_data
from kotte_model import flux_ident_2_data_combination
from process_ident_data import process_info_sample
from plot_ident_results import data_utility_plot
from plot_ident_results import plot_parameter_values
from plot_ident_results import parameter_identifibaility_plot
# from plot_ident_results import parameter_experiment_info_plot
from plot_ident_results import parameter_experiment_info_spider


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

# get experimental system steady state data without noise using Convenience kinetics for v3 (kinetics = 2)
exp_xss, exp_fss, exp_ssid, perturbation_details = \
    generate_expdata(y0, cvode_options, ode_parameter_values, noise=0, kinetics=2, dynamic_plot=0,
                     perturbation_plot=0)

# arrange experimental data to form multiple data sets
exp_flux_index = np.array([0, 3, 2, 4, 1, 5])

# get combination of 2 experiments and perform identifiability on all fluxes that require 2 data sets
# get combinations of experimental datasets
experimental_datasets_2_expts, \
    experiment_choice, combination_choice = arrange_experimental_data(exp_xss, exp_fss, perturbation_details,
                                                                      experiments_per_set=2, flux_id=exp_flux_index,
                                                                      experiment_choice=[0, 1, 2, 3, 4, 5,
                                                                                         11, 12, 13, 14, 15,
                                                                                         16, 17, 18, 19, 20])
# different types of experiments 0 - wt, perturbations: 1 - acetate, 2 - k1cat, 3 - V3max, 4 - V2max
experiment_type_indices = [[0],
                           [1, 2, 3, 4, 5],
                           [6, 7, 8, 9, 10],
                           [11, 12, 13, 14, 15],
                           [16, 17, 18, 19, 20]]

print('Practical Identifiability Analysis of v1 with 2 parameters: k1cat and K1ac\n')
# choose which identifiability functions to test
ident_fun_choice = [0]
# perform identifiability when v1 is written with k1cat*E in the numerator
ident_details_v1_k1cat = flux_ident_2_data_combination(experimental_datasets_2_expts, choose=combination_choice,
                                                       flux_ids=[1], flux_choice=[2], ident_fun_choice=ident_fun_choice)
print('Identifiability analysis of v1 with 2 parameters (k1cat and K1ac) complete.\n')

# data processing - do not combine fluxes

data_list_v1_k1cat, max_parameter_v1_k1cat, true_value_v1_k1cat, experiment_info_v1_k1cat,\
combined_data_list_v1_k1cat, combined_max_parameter_v1_k1cat, combined_true_value_v1_k1cat, \
combined_experiment_info_v1_k1cat = process_info_sample(ident_details_v1_k1cat,
                                                        experimental_datasets_2_expts,
                                                        experiment_type_indices,
                                                        combine_fluxes=0,
                                                        ident_fun_choice=ident_fun_choice)

# plot parameter identifibaility for all fluxes using 2 data combinations
parameter_identifibaility_plot(max_parameter_v1_k1cat)
# plot experiment type in each position based on all parameter
# identifiable data combinations for each parameter
# parameter_experiment_info_plot(experiment_info_2)
parameter_experiment_info_spider(experiment_info_v1_k1cat)
# plot true parameter values and determined parameter values
plot_parameter_values(true_value_v1_k1cat)

# plot utility of data sets (number of data sets identifying n, n-1, n-2, ...., 1, 0 parameters
data_utility_plot(data_list_v1_k1cat)


print("\n Run Complete \n")
