from create_experiment_data import retrieve_experimental_data_from_file
from identifiability_analysis import data_numerical_ident
from numerical_ident import process_opt_solution
from names_strings import true_parameter_values
from validate_estimation import process_validation_info
from plot_ident_results import parameter_values_plot
from plot_ident_results import validation_plot
import os.path


# create data for identifiability analysis
# from create_experiment_data import create_data_for_flux
# create_data_for_flux(flux_id='v3a', noise=0, number_samples=1, kinetics=1)

# extract experimental data from file
new_data_file_name = os.path.join(os.getcwd(), 'exp/exp_v3_4_experiments_mwc')
index_labels = ['sample_name', 'data_set_id', 'experiment_id']
arranged_data_df = retrieve_experimental_data_from_file(data_file_name=new_data_file_name,
                                                        multi_index_label=index_labels)

# perform identifiability when v3 parameters are written using convenience kinetics
# get combination of 3 experiments and perform identifiability on all fluxes that require 3 data sets
storage_file_name = os.path.join(os.getcwd(), 'ident/ident_numerical_v3_mwc_l1_obj')
# lexographic ordering of df indices
all_exp_data = data_numerical_ident(arranged_data_df, 'sample_0')

# NLP solver options
optim_options = {"solver": "ipopt",
                 "opts": {"ipopt.tol": 1e-18}}
# optim_options = {"solver": "sqpmethod",\
#                  "opts": {"qpsol": "qpoases"}}
initial_value = [[80, 80, 80, 80, 0, 0, 0, 0]]  # [100, 400, 400, 0, 0, 0], [50, 50, 400, 0, 0, 0]
# randomized_initial_values = generate_random_initial_conditions(initial_value, 10, negative=1)
problem = {"lbx": 8 * [0],
           "ubx": [2, 1, 1, 5, 1e-20, 1e-20, 1e-20, 1e-20],
           "lbg": 4 * [0],
           "ubg": 4 * [0]}
from numerical_ident import solve_multiple_initial_conditions
all_sol_df, _ = solve_multiple_initial_conditions(all_initial_conditions=initial_value,
                                                  experimental_data=all_exp_data, chosen_fun=2, prob=problem,
                                                  optim_options=optim_options, number_of_parameters=4, flux_id=3,
                                                  flux_choice=[3], exp_df=arranged_data_df,
                                                  file_name=storage_file_name)

index_labels = ['sample_name', 'data_set_id']
numerical_ident_df = retrieve_experimental_data_from_file(data_file_name=storage_file_name,
                                                          multi_index_label=index_labels)
all_parameter_info = process_opt_solution(numerical_ident_df, arranged_data_df, [], [], [], [])

# get default parameter values
default_parameter_values = true_parameter_values()

# extract all parameter values
# from process_ident_data import get_parameter_value
# validation_info = get_parameter_value(all_parameter_info, numerical_ident_df)
# initial value used to generate experimental data
# import numpy as np
# y0 = np.array([5, 1, 1])
# # integrator options
# cvode_options = {'iter': 'Newton', 'discr': 'Adams', 'atol': 1e-10, 'rtol': 1e-10, 'time_points': 200,
#                  'display_progress': False, 'verbosity': 50}
# # validate all parameter values
validation_file_name = os.path.join(os.getcwd(), 'validate/ident_numerical_validate_v3_l1_obj')
# from validate_estimation import validate_model
# validate_model(y0, cvode_options, default_parameter_values, validation_info,
#                save_file_name=validation_file_name,
#                ss=1, dyn=0, noise=1, kinetics=2)

# get parameter value plot
parameter_values_plot(all_parameter_info, default_parameter_values, box=False, violin=True)

# retrieve validation info from file
validate_index_labels = ['estimate_id', 'sample_name', 'data_set_id', 'experiment_id']
validate_df = retrieve_experimental_data_from_file(data_file_name=validation_file_name,
                                                   multi_index_label=validate_index_labels)
# validation plot
original_experiment_file = os.path.join(os.getcwd(), 'exp/experiments')
exp_df = retrieve_experimental_data_from_file(data_file_name=original_experiment_file,
                                              multi_index_label=['sample_name', 'experiment_id'])

all_c_info, all_f_info = process_validation_info(validate_df, exp_df)
validation_plot({"concentration": all_c_info, "flux": all_f_info}, flux=True, flux_id=['v1', 'v2', 'v3', 'v5'],
                experiment_dist=True)

print("Run complete\n")
