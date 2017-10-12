# generate experimental ss/dynamic data by solving model ode using assimulo for different input acetate concentrations
import matplotlib.pyplot as plt
# define rhs function - defined in kotte_model.py
from kotte_model import *
from simulate_ode import run_ode_sims

# main function
if __name__=='__main__':

    y0 = np.array([5, 1, 1])
    cvode_options = ['Newton', 'Adams', 1e-6, 1e-6, 100]
    ode_par_val = np.array([.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, .1])
    cvode_options.append(ode_par_val)

    time, y_dynamic = run_ode_sims(kotte_ode, y0, cvode_options, 100)[:2]
    # calculate and plot dynamic flux data
    flux_function = lambda x: kotte_flux(x, ode_par_val)
    flux_dynamic = map(flux_function, y_dynamic)
    plt.plot(time, flux_dynamic, color="b")
    plt.show()

    time_ck, y_dynamic_ck = run_ode_sims(kotte_ck_ode, y0, cvode_options, 100)[:2]
    # calculate and plot dynamic ck flux data
    flux_function_ck = lambda x: kotte_ck_flux(x, ode_par_val)
    flux_dynamic_ck = map(flux_function_ck, y_dynamic_ck)
    plt.plot(time, flux_dynamic_ck, color="g")
    plt.show()
