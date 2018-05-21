import numpy as np
from create_experiment_data import retrieve_experimental_data_from_file
from process_ident_data import process_ident
from process_ident_data import get_parameter_value
from names_strings import true_parameter_values
from validate_estimation import validate_model
from plot_ident_results import data_utility_plot
from plot_ident_results import plot_parameter_values


# create data for identifiability analysis
# from create_experiment_data import create_data_for_flux
# create_data_for_flux(flux_id='v3', noise=1, number_samples=5)

# extract data from file
new_data_file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels' \
                     '\ident\python2\ss-ident\exp_v3_3_experiments_noise_5_samples'
index_labels = ['sample_name', 'data_set_id', 'experiment_id']
arranged_data_df = retrieve_experimental_data_from_file(data_file_name=new_data_file_name,
                                                        multi_index_label=index_labels)

# perform identifiability when v3 parameters are written for root (1)
# get combination of 3 experiments and perform identifiability on all fluxes that require 3 data sets
# print('Practical Identifiability Analysis of v3 with 3 parameters: V3max, K3fdp and K3pep \n')
storage_file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels' \
                     '\ident\python2\ss-ident\ident_v3_noise_root_1'
# choose identifiability functions to test
ident_fun_choice = [0]
# test identifiability and store data to file
# from kotte_model import flux_ident_3_data_combination
# flux_ident_3_data_combination(arranged_data_df, flux_ids=[3], flux_choice=[1],
#                               ident_fun_choice=ident_fun_choice, file_name=storage_file_name)

# retrieve identifiability data from file and process information
ident_index_label = ['sample_name', 'data_set_id']
# retrieve identifiability info from file
ident_df = retrieve_experimental_data_from_file(storage_file_name, ident_index_label)
all_parameter_info = process_ident(ident_df, arranged_data_df)

# run dynamic simulations to obtain ss data based on estimated parameter values
# get info from data sets that identify all 3 parameters
validation_info = get_parameter_value(all_parameter_info, ident_df)
# get default parameter values
default_parameter_values = true_parameter_values()
# initial value used to generate experimental data
y0 = np.array([5, 1, 1])
# integrator options
cvode_options = ('Newton', 'Adams', 1e-10, 1e-10, 200)
validate_model(y0, cvode_options, default_parameter_values, validation_info,
               save_file_name='C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\ident\python2' \
                              '\ss-ident\ident_validate_v3_noise_root_1',
               ss=1, dyn=0, noise=0, kinetics=2, target_data=range(0, 5))

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
# data_list_v3_root2, max_parameter_v3_root2, true_value_v3_root2, experiment_info_v3_root2, \
#     combined_data_list_v3_root2, combined_max_parameter_v3_root2, combined_true_value_v3_root2, \
#     combined_experiment_info_v3_root2 = process_info_sample(ident_details_v3_root2,
#                                                             experimental_datasets_3_expts,
#                                                             experiment_type_indices,
#                                                             ident_fun_choice=ident_fun_choice)


# plot parameter identifibaility for all fluxes using 2 data combinations
# parameter_identifibaility_plot(max_parameter_v3_root2, noise=1)
# plot experiment type in each position based on all parameter
# identifiable data combinations for each parameter
# parameter_experiment_info_spider(experiment_info_v3_root2, noise=1)
# plot true parameter values and determined parameter values
plot_parameter_values(true_value_v3_root2, noise=1, data_sets_to_plot=range(0, 2000))
# plot utility of data sets (number of data sets identifying n, n-1, n-2, ...., 1, 0 parameters
data_utility_plot(data_list_v3_root2, noise=1)

print("\n Run Complete \n")