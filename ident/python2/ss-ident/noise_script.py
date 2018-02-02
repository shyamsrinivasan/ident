import numpy as np
from generate_expdata import generate_expdata
from simulate_data import arrange_experimental_data
from kotte_model import flux_ident_2_data_combination
from kotte_model import flux_ident_3_data_combination
from process_ident_data import process_info_sample
from plot_ident_results import data_utility_plot
from plot_ident_results import parameter_identifibaility_plot
from plot_ident_results import parameter_experiment_info_plot


# generate noisy experimental data for testing identifiability
y0 = np.array([5, 1, 1])
# default parameter values
cvode_options = ('Newton', 'Adams', 1e-10, 1e-10, 200)
ode_parameter_values = np.array([.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, .1])

# get experimental system steady state data without noise
exp_xss, exp_fss, exp_ssid, perturbation_details = \
    generate_expdata(y0, cvode_options, ode_parameter_values, noise=1, number_of_samples=50)

# arrange experimental data to form multiple data sets
exp_flux_index = np.array([0, 3, 2, 4])

# get combination of 2 experiments and perform identifiability on all fluxes that require 2 data sets
print('Practical Identifiability Analysis of fluxes with 2 parameters \n')
choose_2 = range(0, 306) # choose numbr of experimental datasets to use of analysis
# get combinations of experimental datasets
experimental_datasets_2_expts = \
    arrange_experimental_data(exp_xss, exp_fss, perturbation_details, 2, exp_flux_index, choose_2)
ident_details_2 = flux_ident_2_data_combination(experimental_datasets_2_expts, choose=choose_2,
                                                flux_ids=[1, 2], flux_choice=[1, 0])
print('Identifiability analysis for fluxes with 2 parameters complete.\n')
# data processing
experiment_type_indices = [[0], [1, 2], [3, 4, 5, 6, 7], [8, 9, 10, 11, 12], [13, 14, 15, 16, 17]]
data_list_2, max_parameter_2, experiment_info_2, \
combined_data_list_2, combined_max_parameter_2, \
combined_experiment_info_2 = process_info_sample(ident_details_2,
                                                 experimental_datasets_2_expts,
                                                 experiment_type_indices,
                                                 perturbation_details, combine_fluxes=1)

# plot parameter identifibaility for all fluxes using 2 data combinations
parameter_identifibaility_plot(max_parameter_2, noise=1)
# plot experiment type in each position based on all parameter
# identifiable data combinations for each parameter
parameter_experiment_info_plot(experiment_info_2, noise=1)

# plot utility of data sets (number of data sets identifying n, n-1, n-2, ...., 1, 0 parameters
data_utility_plot(data_list_2, noise=1)

# get combination of 3 experiments and perform identifiability on all fluxes that require 3 data sets
print('Practical Identifiability Analysis of fluxes with 3 parameters \n')
choose_3 = range(0, 4896) # choose numbr of experimental datasets to use of analysis
# get combinations of experimental datasets
experimental_datasets_3_expts = \
    arrange_experimental_data(exp_xss, exp_fss, perturbation_details, 3, exp_flux_index, choose=choose_3)
ident_details_3 = flux_ident_3_data_combination(experimental_datasets_3_expts, choose=choose_3,
                                                flux_ids=[3], flux_choice=[1])
print('Identifiability analysis for fluxes with 2 parameters complete.\n')
# data processing
data_list_3, max_parameter_3, experiment_info_3 = process_info_sample(ident_details_3,
                                                                      experimental_datasets_3_expts,
                                                                      experiment_type_indices,
                                                                      perturbation_details)

# plot parameter identifibaility for all fluxes using 3 data combinations
parameter_identifibaility_plot(max_parameter_3, noise=1)
# plot experiment type in each position based on all parameter
# identifiable data combinations for each parameter
parameter_experiment_info_plot(experiment_info_3, noise=1)

# plot utility of data sets (number of data sets identifying n, n-1, n-2, ...., 1, 0 parameters
data_utility_plot(data_list_3, noise=1)

print("Run complete\n")
