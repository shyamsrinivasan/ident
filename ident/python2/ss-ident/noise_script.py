import numpy as np
from generate_expdata import generate_expdata
from simulate_data import arrange_experimental_data
from kotte_model import establish_kotte_flux_identifiability
from process_ident_data import process_info_sample
from process_ident_data import parameter_plot_data_per_sample
from plot_ident_results import flux_parameter_plot

# generate noisy experimental data for testing identifiability
y0 = np.array([5, 1, 1])
# default parameter values
cvode_options = ('Newton', 'Adams', 1e-10, 1e-10, 200)
ode_parameter_values = np.array([.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, .1])

# get experimental system steady state data without noise
exp_xss, exp_fss, exp_ssid, perturbation_details = \
    generate_expdata(y0, cvode_options, ode_parameter_values, noise=1, number_of_samples=25)

# arrange experimental data to form multiple data sets
exp_flux_index = np.array([0, 3, 2, 4])
# choose numbr of experimental datasets for which identifiability is to be calculated
choose = range(0, 4896)
# get combinations of experimental datasets
experimental_datasets = arrange_experimental_data(exp_xss, exp_fss, perturbation_details, 3, exp_flux_index, choose)

# identifiability for all kotte fluxes
ident_details = establish_kotte_flux_identifiability(experimental_datasets, choose=choose)
print('Perturbation analysis for identifiability complete.\n')

# data processing
data_list, original_data_ident, combo_data_ident, max_parameter = process_info_sample(ident_details,
                                                                                      experimental_datasets,
                                                                                      perturbation_details)
# get data for plots
total_ident_data, fraction_ident_data, all_boolean_p_id = parameter_plot_data_per_sample(original_data_ident, 1)
# plot
flux_parameter_plot(total_ident_data, fraction_ident_data)
# flux_parameter_plot(fraction_ident_data)

# get different classes of datasets (containing different experiments)
experiment_sets = [[0], [1, 2], [3, 4, 5, 6, 7], [8, 9, 10, 11, 12], [13, 14, 15, 16, 17]]
