import numpy as np
from generate_expdata import generate_expdata
from kotte_model import establish_kotte_flux_identifiability
from kotte_model import arrange_experimental_data
from process_ident_data import process_info
from process_ident_data import flux_parameter_plot_data

from plot_ident_results import flux_parameter_plot

# generate noisy experimental data for testing identifiability
y0 = np.array([5, 1, 1])
# default parameter values
cvode_options = ('Newton', 'Adams', 1e-10, 1e-10, 200)
ode_parameter_values = np.array([.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, .1])

# get experimental system steady state data without noise
exp_xss, exp_fss, exp_ssid, perturbation_details = \
    generate_expdata(y0, cvode_options, ode_parameter_values, noise=1, number_of_samples=2)

# arrange experimental data to form multiple data sets
exp_flux_index = np.array([0, 3, 2, 4])
