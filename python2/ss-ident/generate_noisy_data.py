# generate noisy data for different types of kinetics used for regulated metabolic reactions
from kotte_model import *
from simulate_ode import run_ode_sims
from add_noise import add_noise_dynamic


def generate_data(kinetics):
    y0 = np.array([5, 1, 1])
    cvode_options = ['Newton', 'Adams', 1e-6, 1e-6, 100]
    ode_par_val = np.array([.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, .1])
    cvode_options.append(ode_par_val)
    if kinetics == 1:  # MWC kinetics
        time, y_dynamic = run_ode_sims(kotte_ode, y0, cvode_options, 100)[:2]
        # calculate dynamic flux data
        flux_dynamic = np.array(map(lambda x: kotte_flux(x, ode_par_val), y_dynamic))

    elif kinetics == 2:  # Convenience kinetics
        time, y_dynamic = run_ode_sims(kotte_ck_ode, y0, cvode_options, 100)[:2]
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


def generate_noisy_data(kinetics):
    tout, concentration_steady_state, flux_steady_state, \
    concentration_dynamic, flux_dynamic = generate_data(kinetics)

    # add noise to dynamic data
    noisy_concentration_dynamic, noisy_flux_dynamic = add_noise_dynamic(concentration_dynamic, flux_dynamic)
    noisy_concentration_steady_state = noisy_concentration_dynamic[-1, :]
    noisy_flux_steady_state = noisy_flux_dynamic[-1, :]

    return noisy_concentration_steady_state, noisy_flux_steady_state, noisy_concentration_dynamic, noisy_flux_dynamic