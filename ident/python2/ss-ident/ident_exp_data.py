from create_experiment_data import create_data_for_flux


# create data for identifiability analysis (w/o noise)
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
