import numpy as np
import scipy.optimize as opt
from generate_expdata import initialize_to_ss
from kotte_model import kotte_ck_flux
from kotte_model import kotte_ode
import matplotlib.pyplot as plt
from sys import path
# path.append(r"C:\Users\shyam\Anaconda3\envs\py27\Lib\site-packages\casadi-py27-657a052")
# from casadi import *


def kotte_nlae(y, p):
    yres = kotte_ode(0, y, p)
    return yres


def solve_nlae(x0, model_parameters):
    # x_sol, info, flag, err_msg = opt.fsolve(kotte_nlae, x0, model_parameters, xtol=1e-3, full_output=1)
    options = {'ftol': 1e-8, 'disp': 1, 'maxiter': 5000}
    wrapped_fun = lambda y: kotte_nlae(y, model_parameters)
    sol = opt.newton_krylov(wrapped_fun, x0, maxiter=10000, verbose=1, f_tol=1e-8)
    return sol


def parameter_exp(flux, x, ac):
    v1max = flux[0, 0] * flux[0, 1] * (ac[0] - ac[1])/(-(ac[1] * flux[0, 0] - ac[0] * flux[0, 1]))
    # k1ac_no_enzyme_numerator_value = ac1 * (ac2 * v11 - ac2 * v12)
    # k1ac_no_enzyme_denominator_value = -ac2 * v11 + ac1 * v12
    # k1ac
    return v1max


if __name__ == '__main__':
    # initial value of concentrations
    y0 = np.array([5, 1, 1])
    # initial parameter values
    parameter_values = np.array([.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, .1])
    # get initial steady state from ode solution
    cvode_options = ('Newton', 'Adams', 1e-10, 1e-10, 2000)
    initial_ss = initialize_to_ss(y0, cvode_options, parameter_values, noise=0)

    # solve nlae using fsolve strating from initial_ss
    x = solve_nlae(initial_ss["y"], parameter_values)
    f, ax = plt.subplots(1, 1)
    # ax.plot()