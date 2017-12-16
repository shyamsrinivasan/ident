import numpy as np
from generate_expdata import generate_expdata
from kotte_model import establish_kotte_flux_identifiability
from kotte_model import arrange_experimental_data
from process_ident_data import process_info
# from kotte_model import write_results_2_file
# from plot_ident_results import flux_parameter_plot_data
from process_ident_data import useful_experiments


# generate noisy experimental data for testing identifiability
y0 = np.array([5, 1, 1])
# default parameter values
cvode_options = ('Newton', 'Adams', 1e-10, 1e-10, 200)
ode_parameter_values = np.array([.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, .1])

# get noisy experimental system steady state data
noisy_exp_xss, noisy_exp_fss, noisy_exp_ssid, perturbation_details = \
    generate_expdata(y0, cvode_options, ode_parameter_values)

# arrange experimental data to form multiple data sets
exp_flux_index = np.array([0, 3, 2, 4])
# get combinations of experimental datasets
experimental_datasets = arrange_experimental_data(noisy_exp_xss, noisy_exp_fss, perturbation_details, 3, exp_flux_index)

# identifiability for all kotte fluxes
ident_details = establish_kotte_flux_identifiability(experimental_datasets, choose=10)
print('Perturbation analysis for identifiability complete.\n')

# data processing
data_list, original_data_ident, combo_data_ident, max_parameter = process_info(ident_details,
                                                                               experimental_datasets,
                                                                               perturbation_details)

# plot results
file_destination = 'C:\Users\shyam\Documents\Courses\CHE1125Project\Results\ident\python2\\figure_1'
# plot parameters for each flux and the number of data sets that enable their identification
# flux_parameter_plot_data(original_data_ident)
# plot details on experiments in identifiable data sets
all_parameter_exp_id = useful_experiments(original_data_ident)

# save plot
#if file_destination:
#    # save figure to file as png and eps
#    plt.savefig(file_destination + '.eps', format='png', dpi=2000)
#    plt.savefig(file_destination + '.png', format='eps', dpi=2000)


# create data for write_2_file and write to file
# write_results_2_file(ident_details, 3, fp_list, data_list)


# clear workspace (removes all module names and objects)
# import sys
# sys.modules[__name__].__dict__.clear()
