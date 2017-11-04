# generate experimental ss/dynamic data by solving model ode using assimulo for different input acetate concentrations
import matplotlib.pyplot as plt
import numpy as np
from generate_noisy_data import generate_noisy_data
from plot_profiles import plot_dynamic_course


# main function
if __name__ == '__main__':

    y0 = np.array([5, 1, 1])
    # generate data using MWC Kientics
    time, y_noisy_steady_state, flux_noisy_steady_state, \
    y_noisy_dynamic, flux_noisy_dynamic, \
    y_steady_state, flux_steady_state, \
    y_dynamic, flux_dynamic = generate_noisy_data(y0, 1)

    # generate data using Convenience Kinetics
    time_ck, y_ck_noisy_steady_state, flux_ck_noisy_steady_state, \
    y_ck_noisy_dynamic, flux_ck_noisy_dynamic, \
    y_ck_steady_state, flux__ck_steady_state, \
    y_ck_dynamic, flux_ck_dynamic = generate_noisy_data(y0, 2)

    # plot ck dynamic data
    plot_dynamic_course(time_ck, y_ck_dynamic, flux_ck_dynamic, 3)
    plot_dynamic_course(time_ck, y_ck_noisy_dynamic, flux_ck_noisy_dynamic, 3)

    #
    f, axx = plt.subplots(2, 2, sharex='col', sharey='row')
    axx[0, 0].plot(time, y_dynamic, color='r')
    axx[0, 0].set_title('Concentrations, MWC Kinetics')
    axx[0, 1].plot(time, y_ck_dynamic, color='r')
    axx[0, 1].set_title('Concentrations, Convenience Kinetics')
    axx[1, 0].plot(time, flux_dynamic, color='b')
    axx[1, 0].set_title('Flux, MWC Kinetics')
    axx[1, 1].plot(time, flux_ck_dynamic, color='b')
    axx[1, 1].set_title('Flux, Convenience Kinetics')
    f.subplots_adjust(hspace=.3)
    plt.show()
