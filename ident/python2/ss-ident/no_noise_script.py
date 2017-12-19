import numpy as np
from generate_expdata import generate_expdata
from kotte_model import establish_kotte_flux_identifiability
from kotte_model import arrange_experimental_data
from process_ident_data import process_info
from process_ident_data import get_flux_plots_data
from plot_ident_results import flux_parameter_plot

# generate noisy experimental data for testing identifiability
y0 = np.array([5, 1, 1])
# default parameter values
cvode_options = ('Newton', 'Adams', 1e-10, 1e-10, 200)
ode_parameter_values = np.array([.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, .1])

# get experimental system steady state data without noise
exp_xss, exp_fss, exp_ssid, perturbation_details = \
    generate_expdata(y0, cvode_options, ode_parameter_values, noise=0)

# arrange experimental data to form multiple data sets
exp_flux_index = np.array([0, 3, 2, 4])
# get combinations of experimental datasets
experimental_datasets = arrange_experimental_data(exp_xss, exp_fss, perturbation_details, 3, exp_flux_index)

# identifiability for all kotte fluxes
ident_details = establish_kotte_flux_identifiability(experimental_datasets, choose=500)
print('Perturbation analysis for identifiability complete.\n')

# data processing
data_list, original_data_ident, combo_data_ident, max_parameter = process_info(ident_details,
                                                                               experimental_datasets,
                                                                               perturbation_details)

# plot results
# file_destination = 'C:\Users\shyam\Documents\Courses\CHE1125Project\Results\ident\python2\\figure_1'
# plot parameters for each flux and the number of data sets that enable their identification
# get data for plots
total_ident_data, fraction_ident_data, _ = get_flux_plots_data(original_data_ident, 1)
# plot
flux_parameter_plot(total_ident_data)
flux_parameter_plot(fraction_ident_data)

# get different classes of datasets (containing different experiments)
experiment_sets = [[0, 1, 2], [3, 4, 5, 6, 7], [8, 9, 10, 11, 12], [13, 14, 15, 16, 17]]


