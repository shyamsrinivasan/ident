from create_experiment_data import retrieve_experimental_data_from_file
from process_ident_data import process_ident
from validate_estimation import process_validation_info
from plot_ident_results import exp_info_plot
from plot_ident_results import identifiability_plot
from plot_ident_results import parameter_values_plot
from plot_ident_results import validation_plot
from names_strings import true_parameter_values


# create data for identifiability analysis
# from create_experiment_data import create_data_for_flux
# create_data_for_flux(flux_id='v5', noise=0, number_samples=1, kinetics=1)

# extract experimental data from file
new_data_file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels' \
                     '\ident\python2\ss-ident\exp_v5_2_experiments_mwc'
index_labels = ['sample_name', 'data_set_id', 'experiment_id']
arranged_data_df = retrieve_experimental_data_from_file(data_file_name=new_data_file_name,
                                                        multi_index_label=index_labels)

# perform identifiability when v1 parameters with combination of 2 experiments
# print('Practical Identifiability Analysis of v1 with 2 parameters: V1max and K1ac\n')
storage_file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels' \
                     '\ident\python2\ss-ident\ident_v5_root_2_mwc'
# choose identifiability functions to test
# ident_fun_choice = [2]
# # test identifiability and store data to file
# from kotte_model import flux_ident_2_data_combination
# flux_ident_2_data_combination(arranged_data_df, flux_ids=[5], flux_choice=[2],
#                               ident_fun_choice=ident_fun_choice, file_name=storage_file_name)

# retrieve identifiability data from file and process information
ident_index_label = ['sample_name', 'data_set_id']
# retrieve identifiability info from file
ident_df = retrieve_experimental_data_from_file(storage_file_name, ident_index_label)
all_parameter_info = process_ident(ident_df, arranged_data_df)

# get default parameter values
default_parameter_values = true_parameter_values()
# # get parameter value plot
parameter_values_plot(all_parameter_info, default_parameter_values, violin=True, box=False)

# get identifiability plot
identifiability_plot(all_parameter_info)

# get experiment info plot
exp_info_plot(all_parameter_info)

# run dynamic simulations to obtain ss data based on estimated parameter values
# get info from data sets that identify all 3 parameters
validation_file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\ident\python2' \
                       '\ss-ident\ident_validate_v5_root_2_mwc'
# from process_ident_data import get_parameter_value
# validation_info = get_parameter_value(all_parameter_info, ident_df)

# initial value used to generate experimental data
# import numpy as np
# y0 = np.array([5, 1, 1])
# from validate_estimation import validate_model
# # integrator options
# cvode_options = {'iter': 'Newton', 'discr': 'Adams', 'atol': 1e-10, 'rtol': 1e-10, 'time_points': 200,
#                  'display_progress': False, 'verbosity': 50}
# validate_model(y0, cvode_options, default_parameter_values, validation_info,
#                save_file_name=validation_file_name,
#                ss=1, dyn=0, noise=0, kinetics=1)

# retrieve validation info from file
validate_index_labels = ['estimate_id', 'sample_name', 'data_set_id', 'experiment_id']
validate_df = retrieve_experimental_data_from_file(data_file_name=validation_file_name,
                                                   multi_index_label=validate_index_labels)
# validation plot
original_experiment_file = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\ident\python2' \
                            '\ss-ident\experiments_mwc'
exp_df = retrieve_experimental_data_from_file(data_file_name=original_experiment_file,
                                              multi_index_label=['sample_name', 'experiment_id'])

all_c_info, all_f_info = process_validation_info(validate_df, exp_df)
validation_plot({"concentration": all_c_info, "flux": all_f_info}, flux=True, flux_id=['v1', 'v2', 'v3', 'v5'],
                experiment_dist=False)

print("Run complete\n")
