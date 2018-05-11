import numpy as np
from create_experiment_data import retrieve_experimental_data
from simulate_data import arrange_experimental_data
from kotte_model import flux_ident_3_data_combination
from process_ident_data import process_info_sample
from plot_ident_results import parameter_identifibaility_plot
from plot_ident_results import data_utility_plot
from plot_ident_results import plot_parameter_values
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

multi_index_labels = ['sample_name', 'experiment_id']
file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\ident\python2\ss-ident\experiments_noise_5_samples'
exp_ss, perturbation_details = retrieve_experimental_data(file_name=file_name, multi_index_lablel=multi_index_labels)

# arrange experimental data to form multiple data sets
exp_flux_index = np.array([0, 3, 2, 4, 1, 5])

# get combination of 3 experiments and perform identifiability on all fluxes that require 3 data sets
print('Practical Identifiability Analysis of v3 with 3 parameters: V3max, K3fdp and K3pep \n')
# choose identifiability functions to test
ident_fun_choice = [0]
# get combinations of experimental datasets
experimental_datasets_3_expts, \
    experiment_choice, combination_choice = arrange_experimental_data(exp_ss["y"], exp_ss["flux"], perturbation_details,
                                                                      experiments_per_set=3,
                                                                      experiment_choice=[0,
                                                                                         1, 2, 3, 4, 5,
                                                                                         6, 7, 8, 9, 10,
                                                                                         16, 17, 18, 19, 20])
# different types of experiments 0 - wt, perturbations: 1 - acetate, 2 - k1cat, 3 - V3max, 4 - V2max
experiment_type_indices = [[0],
                           [1, 2, 3, 4, 5],
                           [6, 7, 8, 9, 10],
                           [11, 12, 13, 14, 15],
                           [16, 17, 18, 19, 20]]
# perform identifiability when v3 parameters are written for root (1)
ident_details_v3_root1 = flux_ident_3_data_combination(experimental_datasets_3_expts, choose=combination_choice,
                                                       flux_ids=[3], flux_choice=[1], ident_fun_choice=ident_fun_choice)
print('Identifiability analysis for v3 with 3 parameters (V3max, K3fdp and K3pep) complete.\n')

# data processing
data_list_v3_root1, max_parameter_v3_root1, true_value_v3_root1, experiment_info_v3_root1, \
    combined_data_list_v3_root1, combined_max_parameter_v3_root1, combined_true_value_v3_root1, \
    combined_experiment_info_v3_root1 = process_info_sample(ident_details_v3_root1,
                                                            experimental_datasets_3_expts,
                                                            experiment_type_indices,
                                                            ident_fun_choice=ident_fun_choice,
                                                            combine_fluxes=0)

# plot parameter identifibaility for all fluxes using 2 data combinations
parameter_identifibaility_plot(max_parameter_v3_root1, noise=1)
# plot experiment type in each position based on all parameter
# identifiable data combinations for each parameter
parameter_experiment_info_spider(experiment_info_v3_root1, noise=1)
from process_ident_data import extract_parameter_values
extract_parameter_values(true_value_v3_root1)
# plot true parameter values and determined parameter values
plot_parameter_values(true_value_v3_root1, noise=1, data_sets_to_plot=range(0, 300))
plot_parameter_values(true_value_v3_root1, noise=1, data_sets_to_plot=range(300, 600))
plot_parameter_values(true_value_v3_root1, noise=1, data_sets_to_plot=range(600, 1000))
plot_parameter_values(true_value_v3_root1, noise=1, data_sets_to_plot=range(1000, 1300))
plot_parameter_values(true_value_v3_root1, noise=1, data_sets_to_plot=range(1300, 1600))
plot_parameter_values(true_value_v3_root1, noise=1, data_sets_to_plot=range(1600, 2000))
plot_parameter_values(true_value_v3_root1, noise=1, data_sets_to_plot=range(2000, 2300))
plot_parameter_values(true_value_v3_root1, noise=1, data_sets_to_plot=range(2300, 2500))
# plot utility of data sets (number of data sets identifying n, n-1, n-2, ...., 1, 0 parameters
data_utility_plot(data_list_v3_root1, noise=1)

print('Practical Identifiability Analysis of v3 with 3 parameters: V3max, K3fdp and K3pep \n')
# perform identifiability when v3 parameters are written for root (2)
ident_details_v3_root2 = flux_ident_3_data_combination(experimental_datasets_3_expts, choose=combination_choice,
                                                       flux_ids=[3], flux_choice=[2], ident_fun_choice=ident_fun_choice)
print('Identifiability analysis for v3 with 3 parameters (V3max, K3fdp and K3pep) complete.\n')
# data processing
data_list_v3_root2, max_parameter_v3_root2, true_value_v3_root2, experiment_info_v3_root2, \
    combined_data_list_v3_root2, combined_max_parameter_v3_root2, combined_true_value_v3_root2, \
    combined_experiment_info_v3_root2 = process_info_sample(ident_details_v3_root2,
                                                            experimental_datasets_3_expts,
                                                            experiment_type_indices,
                                                            ident_fun_choice=ident_fun_choice)


# plot parameter identifibaility for all fluxes using 2 data combinations
parameter_identifibaility_plot(max_parameter_v3_root2, noise=1)
# plot experiment type in each position based on all parameter
# identifiable data combinations for each parameter
parameter_experiment_info_spider(experiment_info_v3_root2, noise=1)
# plot true parameter values and determined parameter values
plot_parameter_values(true_value_v3_root2, noise=1, data_sets_to_plot=range(0, 2000))
# plot utility of data sets (number of data sets identifying n, n-1, n-2, ...., 1, 0 parameters
data_utility_plot(data_list_v3_root2, noise=1)

print("\n Run Complete \n")