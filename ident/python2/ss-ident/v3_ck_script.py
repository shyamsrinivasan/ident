import numpy as np
from create_experiment_data import retrieve_experimental_data_from_file
from kotte_model import flux_ident_3_data_combination
from process_ident_data import process_ident
from process_ident_data import process_info_sample
from plot_ident_results import exp_info_plot
from plot_ident_results import data_utility_plot
from plot_ident_results import plot_parameter_values
from plot_ident_results import plot_parameter_value_hist
from plot_ident_results import parameter_identifibaility_plot


# create data for identifiability analysis
# from create_experiment_data import create_data_for_flux
# create_data_for_flux(flux_id='v3', noise=0, number_samples=1)

# extract experimental data from file
new_data_file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels' \
                     '\ident\python2\ss-ident\exp_v3_3_experiments'
index_labels = ['sample_name', 'data_set_id', 'experiment_id']
arranged_data_df = retrieve_experimental_data_from_file(data_file_name=new_data_file_name,
                                                        multi_index_label=index_labels)

# perform identifiability when v3 parameters are written for root (1)
# get combination of 3 experiments and perform identifiability on all fluxes that require 3 data sets
# print('Practical Identifiability Analysis of v3 with 3 parameters: V3max, K3fdp and K3pep \n')
storage_file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels' \
                     '\ident\python2\ss-ident\ident_v3_root_1'
# choose identifiability functions to test
ident_fun_choice = [0]
# test identifiability and store data to file
# flux_ident_3_data_combination(arranged_data_df, arranged_data_df, flux_ids=[3], flux_choice=[1],
#                               ident_fun_choice=ident_fun_choice, file_name=storage_file_name)

# retrieve identifiability data from file and process information
ident_index_label = ['sample_name', 'data_set_id']
# retrieve identifiability info from file
ident_df = retrieve_experimental_data_from_file(storage_file_name, ident_index_label)
all_parameter_info = process_ident(ident_df, arranged_data_df)

# get identifiability plot
# get experiment info plot
exp_info_plot(all_parameter_info)
# get parameter value plot
# validate model

# different types of experiments 0 - wt, perturbations: 1 - acetate, 2 - k1cat, 3 - V3max, 4 - V2max
experiment_type_indices = [[0],
                           [1, 2, 3, 4, 5],
                           [6, 7, 8, 9, 10],
                           [11, 12, 13, 14, 15],
                           [16, 17, 18, 19, 20]]

# data processing
data_list_v3_root1, max_parameter_v3_root1, true_value_v3_root1, experiment_info_v3_root1, \
    combined_data_list_v3_root1, combined_max_parameter_v3_root1, combined_true_value_v3_root1, \
    combined_experiment_info_v3_root1 = process_info_sample(ident_details_v3_root1,
                                                            experimental_datasets_3_expts,
                                                            experiment_type_indices,
                                                            ident_fun_choice=ident_fun_choice)

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



# plot parameter identifibaility for all fluxes using 3 data combinations
parameter_identifibaility_plot(max_parameter_v3_root1)
# plot true parameter values and determined parameter values
plot_parameter_values(true_value_v3_root1)
plot_parameter_value_hist(true_value_v3_root1)
# plot utility of data sets (number of data sets identifying n, n-1, n-2, ...., 1, 0 parameters
data_utility_plot(data_list_v3_root1)
from process_ident_data import extract_parameter_values
parameter_value_info = extract_parameter_values(true_value_v3_root1)
from validate_estimation import validate_model
validate_model(y0, cvode_options, ode_parameter_values, parameter_value_info, exp_info,
               ss=1, dyn=0, noise=0, kinetics=2, target_data=range(0, 10))

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

# plot parameter identifibaility for all fluxes using 3 data combinations
parameter_identifibaility_plot(max_parameter_v3_root2)
# plot true parameter values and determined parameter values
plot_parameter_values(true_value_v3_root2)

# plot utility of data sets (number of data sets identifying n, n-1, n-2, ...., 1, 0 parameters
data_utility_plot(data_list_v3_root2)

print("Run complete\n")
