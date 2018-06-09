from create_experiment_data import retrieve_experimental_data_from_file
from process_ident_data import process_ident
from validate_estimation import process_validation_info
from plot_ident_results import exp_info_plot
from plot_ident_results import identifiability_plot
from plot_ident_results import parameter_values_plot
from plot_ident_results import validation_plot
from names_strings import true_parameter_values
import os.path


# create data for identifiability analysis
# from create_experiment_data import create_data_for_flux
# create_data_for_flux(flux_id='v3b', noise=0, number_samples=1)

# extract experimental data from file
new_data_file_name = os.path.join(os.getcwd(), 'exp/exp_v3_2_experiments')
index_labels = ['sample_name', 'data_set_id', 'experiment_id']
arranged_data_df = retrieve_experimental_data_from_file(data_file_name=new_data_file_name,
                                                        multi_index_label=index_labels)

# perform identifiability when v3 parameters are written for root (1)
# get combination of 2 experiments and perform identifiability on all fluxes that require 2 data sets
# print('Practical Identifiability Analysis of v3 with 2 parameters: V3max and K3pep \n')
storage_file_name = os.path.join(os.getcwd(), 'ident/ident_v3_k3pep')
# choose identifiability functions to test
# ident_fun_choice = [3]
# # test identifiability and store data to file
# from kotte_model import flux_ident_2_data_combination
# flux_ident_2_data_combination(arranged_data_df, flux_ids=[6], flux_choice=[2],
#                               ident_fun_choice=ident_fun_choice, file_name=storage_file_name)

# retrieve identifiability data from file and process information
ident_index_label = ['sample_name', 'data_set_id']
# retrieve identifiability info from file
ident_df = retrieve_experimental_data_from_file(storage_file_name, ident_index_label)
all_parameter_info = process_ident(ident_df, arranged_data_df)

# run dynamic simulations to obtain ss data based on estimated parameter values
# get info from data sets that identify all 3 parameters
validation_file_name = os.path.join(os.getcwd(), 'validate/ident_validate_v3_k3pep')
# from process_ident_data import get_parameter_value
# validation_info = get_parameter_value(all_parameter_info, ident_df)
# get default parameter values
default_parameter_values = true_parameter_values()

# initial value used to generate experimental data
# from validate_estimation import validate_model
# import numpy as np
# y0 = np.array([5, 1, 1])
# # integrator options
# cvode_options = {'iter': 'Newton', 'discr': 'Adams', 'atol': 1e-10, 'rtol': 1e-10, 'time_points': 500,
#                  'display_progress': False, 'verbosity': 50}
# validate_model(y0, cvode_options, default_parameter_values, validation_info,
#                save_file_name=validation_file_name,
#                ss=1, dyn=0, noise=0, kinetics=2)

# retrieve validation info from file
validate_index_labels = ['estimate_id', 'sample_name', 'data_set_id', 'experiment_id']
validate_df = retrieve_experimental_data_from_file(data_file_name=validation_file_name,
                                                   multi_index_label=validate_index_labels)
# validation plot
original_experiment_file = os.path.join(os.getcwd(), 'exp/experiments')
exp_df = retrieve_experimental_data_from_file(data_file_name=original_experiment_file,
                                              multi_index_label=['sample_name', 'experiment_id'])

all_c_info, all_f_info = process_validation_info(validate_df, exp_df)
validation_plot({"concentration": all_c_info, "flux": all_f_info}, flux=True, flux_id=['v1', 'v2', 'v3', 'v5'])

# get parameter value plot
parameter_values_plot(all_parameter_info, default_parameter_values, violin=True, box=False)

# get identifiability plot
identifiability_plot(all_parameter_info)

# get experiment info plot
exp_info_plot(all_parameter_info)

print("Run complete\n")