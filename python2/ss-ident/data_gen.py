# generate experimental ss/dynamic data by solving model ode using assimulo for different input acetate concentrations
import numpy as np
import matplotlib.pyplot as plt
# define rhs function - defined in kotte_model.py
from kotte_model import kotte_ode
from kotte_model import kotte_flux
from simulate_ode import simulate_ode


def run_example(fun, y_initial, opts, t_final=500, args_1=True):
    """run kotte model ode using cvode from assimulo"""

    # def_par_val = np.array([.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, .1])
    time_points, y_dynamic, prob, solver = simulate_ode(fun, y_initial, t_final, opts)
    if args_1:
        plt.plot(time_points, y_dynamic, color="r")
        plt.xlabel('Time')
        plt.ylabel('Dependent Variables')
        plt.show()

    return time_points, y_dynamic, prob, solver

# main function
if __name__=='__main__':

    y0 = np.array([5, 1, 1])
    cvode_options = ['Newton', 'Adams', 1e-6, 1e-6, 100]
    ode_par_val = np.array([.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, 1])
    cvode_options.append(ode_par_val)

    time, y_dynamic = run_example(kotte_ode, y0, cvode_options, 100)[:2]

    # get flux data for all dynamic concentrations in y_dynamic
    # par_val = np.array([.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, .1])

    # plot dynamic flux data
    flux_function = lambda x: kotte_flux(x, ode_par_val)
    flux_dynamic = map(flux_function, y_dynamic)
    plt.plot(time, flux_dynamic, color="b")
    plt.show()