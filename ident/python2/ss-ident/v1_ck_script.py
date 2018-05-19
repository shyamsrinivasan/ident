import numpy as np
from create_experiment_data import retrieve_experimental_data_from_file
from process_ident_data import process_ident
from process_ident_data import get_parameter_value
from validate_estimation import validate_model
from plot_ident_results import exp_info_plot
from plot_ident_results import identifiability_plot
from plot_ident_results import parameter_values_plot
from names_strings import true_parameter_values
from plot_ident_results import data_utility_plot


# create data for identifiability analysis
# from create_experiment_data import create_data_for_flux
# create_data_for_flux(flux_id='v1', noise=0, number_samples=1)

# extract experimental data from file
new_data_file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels' \
                     '\ident\python2\ss-ident\exp_v1_2_experiments'
index_labels = ['sample_name', 'data_set_id', 'experiment_id']
arranged_data_df = retrieve_experimental_data_from_file(data_file_name=new_data_file_name,
                                                        multi_index_label=index_labels)

# perform identifiability when v1 parameters are written with k1cat
# print('Practical Identifiability Analysis of v1 with 2 parameters: k1cat and K1ac\n')
storage_file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels' \
                     '\ident\python2\ss-ident\ident_v1_k1cat'
# choose which identifiability functions to test
ident_fun_choice = [0]
# test identifiability and store data to file
# from kotte_model import flux_ident_2_data_combination
# flux_ident_2_data_combination(arranged_data_df, flux_ids=[1], flux_choice=[2], ident_fun_choice=ident_fun_choice,
#                               file_name=storage_file_name)

# retrieve identifiability data from file and process information
ident_index_label = ['sample_name', 'data_set_id']
# retrieve identifiability info from file
ident_df = retrieve_experimental_data_from_file(storage_file_name, ident_index_label)
all_parameter_info = process_ident(ident_df, arranged_data_df)

# run dynamic simulations to obtain ss data based on estimated parameter values
validation_info = get_parameter_value(all_parameter_info, ident_df)
# get default parameter values
default_parameter_values = true_parameter_values()
# initial value used to generate experimental data
y0 = np.array([5, 1, 1])
# integrator options
cvode_options = ('Newton', 'Adams', 1e-10, 1e-10, 200)
validate_model(y0, cvode_options, default_parameter_values, validation_info,
               save_file_name='C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\ident\python2'
                              '\ss-ident\ident_validate_v1_k1cat',
               ss=1, dyn=0, noise=0, kinetics=2)

# get parameter value plot
parameter_values_plot(all_parameter_info, default_parameter_values)

# get identifiability plot
identifiability_plot(all_parameter_info)

# get experiment info plot
exp_info_plot(all_parameter_info)

# plot utility of data sets (number of data sets identifying n, n-1, n-2, ...., 1, 0 parameters
data_utility_plot(data_list_v1_k1cat)

print("\n Run Complete \n")
