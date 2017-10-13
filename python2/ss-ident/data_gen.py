# generate experimental ss/dynamic data by solving model ode using assimulo for different input acetate concentrations
import matplotlib.pyplot as plt
# define rhs function - defined in kotte_model.py
from kotte_model import *
from simulate_ode import run_ode_sims
from add_noise import add_noise_dynamic
from plot_profiles import plot_dynamic_course

# main function
if __name__=='__main__':

    y0 = np.array([5, 1, 1])
    cvode_options = ['Newton', 'Adams', 1e-6, 1e-6, 100]
    ode_par_val = np.array([.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, .1])
    cvode_options.append(ode_par_val)

    time, y_dynamic = run_ode_sims(kotte_ode, y0, cvode_options, 100)[:2]
    # calculate dynamic flux data
    flux_dynamic = np.array(map(lambda x: kotte_flux(x, ode_par_val), y_dynamic))
    # add noise to dynamic data
    noisy_y_dynamic, noisy_flux_dynamic = add_noise_dynamic(y_dynamic, flux_dynamic)
    # plot no noise data
    plot_dynamic_course(time, y_dynamic, flux_dynamic, 3)
    plot_dynamic_course(time, noisy_y_dynamic, noisy_flux_dynamic, 3)

    # get ss info from dynamic data
    y_steady_state = y_dynamic[-1,:]
    flux_steady_state = flux_dynamic[-1, :]
    # get noisy ss from dynamic data
    y_noisy_steady_state = noisy_y_dynamic[-1, :]
    flux_noisy_steady_state = noisy_flux_dynamic[-1, :]

    time_ck, y_ck_dynamic = run_ode_sims(kotte_ck_ode, y0, cvode_options, 100)[:2]
    # calculate dynamic flux
    flux_ck_dynamic = np.array(map(lambda x: kotte_ck_flux(x, ode_par_val), y_ck_dynamic))
    # add noise to dynamic data
    noisy_y_ck_dynamic, noisy_flux_ck_dynamic = add_noise_dynamic(y_ck_dynamic, flux_ck_dynamic)
    # plot ck dynamic data
    plot_dynamic_course(time_ck, y_ck_dynamic, flux_ck_dynamic, 3)
    plot_dynamic_course(time_ck, noisy_y_ck_dynamic, noisy_flux_ck_dynamic, 3)

    # get ss info from dynamic data
    y_ck_steady_state = y_ck_dynamic[-1, :]
    flux_ck_steady_state = flux_ck_dynamic[-1, :]
    # get noisy ss from dynamic data
    y_ck_noisy_steady_state = noisy_y_ck_dynamic[-1, :]
    flux_ck_noisy_steady_state = noisy_flux_ck_dynamic[-1, :]

    #
    f, axx = plt.subplots(2, 2, sharex='col', sharey='row')
    axx[0,0].plot(time, y_dynamic, color='r')
    axx[0,0].set_title('Concentrations, No Noise')
    axx[0,1].plot(time, noisy_y_dynamic, color='r')
    axx[0,1].set_title('Concentrations, Noisy')
    axx[1,0].plot(time, flux_dynamic, color='b')
    axx[1,0].set_title('Flux, No Noise')
    axx[1,1].plot(time, noisy_flux_dynamic, color='b')
    axx[1,1].set_title('Flux, Noisy')
    f.subplots_adjust(hspace=.3)
    plt.show()

