import numpy as np
from create_experiment_data import retrieve_experimental_data_from_file
from identifiability_analysis import collect_data
from numerical_ident import identify_all_data_sets


# create data for identifiability analysis
# from create_experiment_data import create_data_for_flux
# create_data_for_flux(flux_id='v3', noise=0, number_samples=1)

# extract experimental data from file
new_data_file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels' \
                     '\ident\python2\ss-ident\exp_v3_3_experiments'
index_labels = ['sample_name', 'data_set_id', 'experiment_id']
arranged_data_df = retrieve_experimental_data_from_file(data_file_name=new_data_file_name,
                                                        multi_index_label=index_labels)

# perform identifiability when v3 parameters are written using mwc kinetics
# get combination of 4 experiments and perform identifiability on all fluxes that require 4 data sets
storage_file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels' \
                     '\ident\python2\ss-ident\ident_v3_root_1'
# lexographic ordering of df indices
arranged_data_df.sort_index(level='data_set_id', inplace=True)
arranged_data_df.sort_index(level='sample_name', inplace=True)
all_exp_data = collect_data(arranged_data_df, 'sample_0')

# NLP solver options
optim_options = {"solver": "ipopt",
                 "opts": {"ipopt.tol": 1e-12}}
# optim_options = {"solver": "sqpmethod",
#                  "opts": {"qpsol": "qpoases"}}
initial_value = [.1, .1, .1, 1, 0, 0, 0, 0]
identify_all_data_sets(experimental_data=all_exp_data, chosen_fun=1,
                       optim_options=optim_options, x0=initial_value)

# generate noisy experimental data for testing identifiability
y0 = np.array([5, 1, 1])
# default parameter values
cvode_options = ('Newton', 'Adams', 1e-10, 1e-10, 200)
ode_parameter_values = {"K1ac": np.array([.1]),
                        "K3fdp": np.array([.1]),
                        "L3fdp": np.array([4e6]),
                        "K3pep": np.array([.1]),
                        "K2pep": np.array([.3]),
                        "vemax": np.array([1.1]),
                        "Kefdp": np.array([.45]),
                        "ne": np.array([2]),
                        "d": np.array([.25]),
                        "V4max": np.array([.2]),
                        "k1cat": np.array([1]),
                        "V3max": np.array([1]),
                        "V2max": np.array([1]),
                        "ac": np.array([.1])}

