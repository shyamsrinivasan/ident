import numpy as np
from create_experiment_data import retrieve_experimental_data_from_file
from identifiability_analysis import data_numerical_ident
from names_strings import true_parameter_values
from validate_estimation import validate_model
from plot_ident_results import plot_numerical_parameter_estimates


# create data for identifiability analysis
# from create_experiment_data import create_data_for_flux
# create_data_for_flux(flux_id='v3', noise=0, number_samples=1)

# extract experimental data from file
new_data_file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels' \
                     '\ident\python2\ss-ident\exp_v3_3_experiments'
index_labels = ['sample_name', 'data_set_id', 'experiment_id']
arranged_data_df = retrieve_experimental_data_from_file(data_file_name=new_data_file_name,
                                                        multi_index_label=index_labels)

# perform identifiability when v3 parameters are written using convenience kinetics
# get combination of 3 experiments and perform identifiability on all fluxes that require 3 data sets
storage_file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels' \
                     '\ident\python2\ss-ident\ident_numerical_v3'
# lexographic ordering of df indices
all_exp_data = data_numerical_ident(arranged_data_df, 'sample_0')

# NLP solver options
optim_options = {"solver": "ipopt",
                 "opts": {"ipopt.tol": 1e-12}}
# optim_options = {"solver": "sqpmethod",\
#                  "opts": {"qpsol": "qpoases"}}
initial_value = [80, 80, 400, 0, 0, 0]
# randomized_initial_values = generate_random_initial_conditions(initial_value, 10, negative=1)
from numerical_ident import solve_multiple_initial_conditions
all_sol_df = solve_multiple_initial_conditions(all_initial_conditions=[initial_value],
                                               experimental_data=all_exp_data, chosen_fun=0,
                                               optim_options=optim_options, number_of_parameters=3, flux_id=3,
                                               flux_choice=[3], exp_df=arranged_data_df, file_name=storage_file_name)

index_labels = ['sample_name', 'data_set_id']
numerical_ident_df = retrieve_experimental_data_from_file(data_file_name=storage_file_name,
                                                          multi_index_label=index_labels)

# plot_numerical_parameter_estimates(v3_all_x0_parameter_info[0])
# extract all parameter values
from process_ident_data import get_parameter_value
parameter_value_info = get_parameter_value(numerical_ident_df)
# validate all parameter values
y0 = np.array([5, 1, 1])
# default parameter values
default_parameter_values = true_parameter_values()
cvode_options = ('Newton', 'Adams', 1e-10, 1e-10, 200)
ode_parameter_values = true_parameter_values()
# validate_model(y0, cvode_options, ode_parameter_values, parameter_value_info, exp_info,
#                ss=1, dyn=0, noise=0, kinetics=2, target_data=range(0, 3300))
print("Run complete\n")
