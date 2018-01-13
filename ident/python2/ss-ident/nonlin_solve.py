import numpy as np
import scipy.optimize as opt
from generate_expdata import initialize_to_ss
from kotte_model import kotte_ck_flux
from kotte_model import kotte_ode
import matplotlib.pyplot as plt
from sys import path
path.append(r"C:\Users\shyam\Anaconda3\envs\py27\Lib\site-packages\casadi-py27-657a052")
from casadi import *


def kotte_casadi_nlae(x0, p):
    y0 = SX.sym('y0')
    y1 = SX.sym('y1')
    y2 = SX.sym('y2')

    K1ac, K3fdp, L3fdp, K3pep, K2pep, vemax, Kefdp, ne, d, V4max, k1cat, V3max, V2max, ac = p

    flux_1 = k1cat * y2 * ac / (ac + K1ac)
    flux_2 = vemax * (1 - 1 / (1 + (Kefdp / y1) ** ne))
    fdp_sat = 1 + y1 / K3fdp
    pep_sat = 1 + y0 / K3pep
    flux_3 = V3max * (fdp_sat - 1) * (fdp_sat ** 3) / (fdp_sat ** 4 + L3fdp * (pep_sat ** (-4)))
    flux_4 = V2max * y0 / (y0 + K2pep)
    flux_5 = V4max * y0
    flux_6 = d * y2
    flux = horzcat(flux_1, flux_2, flux_3, flux_4, flux_5, flux_6)
    # flux = np.hstack((flux_1, flux_2, flux_3, flux_4, flux_5, flux_6))

    yres0 = flux[0] - flux[3] - flux[4]
    yres1 = flux[3] - flux[2]
    yres2 = flux[1] - flux[5]
    # yres_exp = horzcat(yres_pep, yres_fdp, yres_e)

    yres = Function('yres', [y0, y1, y2], [yres0, yres1, yres2])
    g = rootfinder('g', 'newton', yres)

    result = g(y0)

    return None


def kotte_nlae(y, p):
    yres = kotte_ode(0, y, p)
    return yres


def solve_nlae(x0, model_parameters):
    # x_sol, info, flag, err_msg = opt.fsolve(kotte_nlae, x0, model_parameters, xtol=1e-3, full_output=1)
    options = {'ftol': 1e-8, 'disp': 1, 'maxiter': 5000}
    wrapped_fun = lambda y: kotte_nlae(y, model_parameters)
    sol = opt.anderson(wrapped_fun, x0, maxiter=100000, f_tol=1e-8)
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
    fixed_parameter_value = np.array([.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, .1])
    # get initial steady state from ode solution
    cvode_options = ('Newton', 'Adams', 1e-10, 1e-10, 2000)
    initial_ss = initialize_to_ss(y0, cvode_options, fixed_parameter_value, noise=0)

    kotte_casadi_nlae(initial_ss["y"], fixed_parameter_value)

    acetate_values = np.arange(.1, 2, 0.01)
    all_solutions = []
    old_ss = initial_ss["y"]
    for iter_id, ivalue in enumerate(acetate_values):
        print("iteration : {}".format(iter_id+1))
        new_parameter_value = fixed_parameter_value[:]
        new_parameter_value[-1] = ivalue
        # solve nlae using fsolve strating from initial_ss
        x = solve_nlae(old_ss, new_parameter_value)
        ols_ss = x
        all_solutions.append(x)

    f, ax = plt.subplots(1, 1)
    # ax.plot(acetate_values, )