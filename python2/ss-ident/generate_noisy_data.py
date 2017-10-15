# generate noisy data for different types of kinetics used for regulated metabolic reactions
from kotte_model import *
from simulate_ode import run_ode_sims
from add_noise import add_noise_dynamic


def generate_data(y0, all_options, kinetics):
    cvode_options, ode_par_val = all_options
    if kinetics == 1:  # MWC kinetics
        time, y_dynamic = run_ode_sims(kotte_ode, y0, all_options, 100)[:2]
        # calculate dynamic flux data
        flux_dynamic = np.array(map(lambda x: kotte_flux(x, ode_par_val), y_dynamic))

    elif kinetics == 2:  # Convenience kinetics
        time, y_dynamic = run_ode_sims(kotte_ck_ode, y0, all_options, 100)[:2]
        # calculate dynamic flux
        flux_dynamic = np.array(map(lambda x: kotte_ck_flux(x, ode_par_val), y_dynamic))
    else:
        time = []
        y_dynamic = []
        flux_dynamic = []

    # get ss info from dynamic data
    y_steady_state = y_dynamic[-1, :]
    flux_steady_state = flux_dynamic[-1, :]

    return time, y_steady_state, flux_steady_state, y_dynamic, flux_dynamic


def generate_noisy_data(y0, all_options, kinetics):
    tout, concentration_steady_state, flux_steady_state, \
    concentration_dynamic, flux_dynamic = generate_data(y0, all_options, kinetics)

    # add noise to dynamic data
    noisy_concentration_dynamic, noisy_flux_dynamic = add_noise_dynamic(concentration_dynamic, flux_dynamic)
    noisy_concentration_steady_state = noisy_concentration_dynamic[-1, :]
    noisy_flux_steady_state = noisy_flux_dynamic[-1, :]

    return tout, noisy_concentration_steady_state, noisy_flux_steady_state, \
           noisy_concentration_dynamic, noisy_flux_dynamic, \
           concentration_steady_state, flux_steady_state, \
           concentration_dynamic, flux_dynamic