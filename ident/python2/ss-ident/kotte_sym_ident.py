import numpy as np
import matplotlib.pyplot as plt
from generate_noisy_data import generate_noisy_data
from generate_noisy_data import run_noisy_parameter_perturbation
from kotte_model import establish_kotte_flux_identifiability
from kotte_model import arrange_experimental_data
from kotte_model import ident_parameter_name
from kotte_model import kotte_parameter_name
from process_ident_data import process_info
from kotte_model import write_results_2_file
from plot_ident_results import flux_parameter_plot_data

# generate noisy experimental data for testing identifiability
y0 = np.array([5, 1, 1])
all_options_exp_1 = []
all_options_exp_2 = []
all_options_exp_3 = []
# default parameter values
cvode_options = ('Newton', 'Adams', 1e-10, 1e-10, 200)
ode_parameter_values = np.array([.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, .1])

# get initial noisy system steady state
initial_options = (cvode_options, ode_parameter_values)
noisy_initial_ss, _, _, _ = generate_noisy_data(y0, initial_options, 1)

# all parameter perturbations
parameter_perturbation = [(14, 0), (14, 4), (14, 9),
                          (11, .1), (11, .5), (11, 1), (11, -.1), (11, -.5),
                          (12, .1), (12, .5), (12, 1), (12, -.1), (12, -.5),
                          (13, .1), (13, .5), (13, 1), (13, -.1), (13, -.5)]
perturbation_options = {'ode_parameters':ode_parameter_values, 'cvode_options':cvode_options}
noisy_ss, noisy_dynamic, perturbation_details, _, dynamic_info = \
    run_noisy_parameter_perturbation(parameter_perturbation, noisy_initial_ss["y"], perturbation_options)
# plot all dynamic courses
# plot_multiple_dynamics(noisy_dynamic)
# plt.close("all")
# plot_multiple_dynamics(dynamic_info)
# plt.close("all")

noisy_exp_xss = []
noisy_exp_fss = []
noisy_exp_ssid = []
for ss_values in noisy_ss:
    noisy_exp_xss.append(ss_values["y"])
    noisy_exp_fss.append(ss_values["flux"])
    noisy_exp_ssid.append(ss_values["ssid"])

# experimental data based on order of inputs for lambdify expressions
exp_flux_index = np.array([0, 3, 2, 4])
# get combinations of experimental datasets
experimental_datasets = arrange_experimental_data(noisy_exp_xss, noisy_exp_fss, perturbation_details, 3, exp_flux_index)

# identifiability for all kotte fluxes
ident_details = establish_kotte_flux_identifiability(experimental_datasets, choose=10)
print('Perturbation analysis for identifiability complete.\n')

# data processing
data_list, original_data_ident, combo_data_ident, max_parameter = \
    process_info(ident_details,
                 experimental_datasets,
                 perturbation_details,
                 ident_parameter_name, kotte_parameter_name)

# plot results
file_destination = 'C:\Users\shyam\Documents\Courses\CHE1125Project\Results\ident\python2\\figure_1'
flux_parameter_plot_data(original_data_ident)
# save plot
if file_destination:
    # save figure to file as png and eps
    plt.savefig(file_destination + '.eps', format='png', dpi=2000)
    plt.savefig(file_destination + '.png', format='eps', dpi=2000)
# plot_identifiable_parameter(max_parameter)

# create data for write_2_file and write to file
write_results_2_file(ident_details, 3, fp_list, data_list)


# clear workspace (removes all module names and objects)
# import sys
# sys.modules[__name__].__dict__.clear()
