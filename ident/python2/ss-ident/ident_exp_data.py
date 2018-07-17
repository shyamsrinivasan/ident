from create_experiment_data import create_experiment_data
from create_experiment_data import create_data_for_flux
import os.path


file_name = os.path.join(os.getcwd(), 'exp/experiments')
experiment_info_df = create_experiment_data(file_name, noise=0, kinetics=2)

file_name = os.path.join(os.getcwd(), 'exp/experiments_mwc')
experiment_info_df = create_experiment_data(file_name, noise=0, kinetics=1)

# file_name = os.path.join(os.getcwd(), 'exp/experiments_noise_5_samples')
# experiment_info_df = create_experiment_data(file_name, noise=1, kinetics=2, number_samples=5, noise_std=0.05)

# file_name = os.path.join(os.getcwd(), 'exp/experiments_noise_500_samples')
# experiment_info_df = create_experiment_data(file_name, noise=1, kinetics=2, number_samples=500, noise_std=0.05)

# create data for identifiability analysis
# v1 - kcat and vmax with ck model
create_data_for_flux(flux_id='v1', noise=0, number_samples=1)

# v1 - kcat and vmax with mwc model
create_data_for_flux(flux_id='v1', noise=0, number_samples=1, kinetics=1)

# v2 - ck model
create_data_for_flux(flux_id='v2', noise=0, number_samples=1)

# v2 - mwc model
create_data_for_flux(flux_id='v2', noise=0, number_samples=1, kinetics=1)

# v3 - ck model
create_data_for_flux(flux_id='v3', noise=0, number_samples=1)

# v3 - mwc model
create_data_for_flux(flux_id='v3', noise=0, number_samples=1, kinetics=1)

# v3 - ck model with 2 experiments for kfdp and kpep
create_data_for_flux(flux_id='v3b', noise=0, number_samples=1)

# v5 - ck model
create_data_for_flux(flux_id='v5', noise=0, number_samples=1)

# v5 - mwc model
create_data_for_flux(flux_id='v5', noise=0, number_samples=1, kinetics=1)
