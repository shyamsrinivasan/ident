import numpy as np
import scipy.optimize as opt
from kotte_model import kotte_ck_flux
import matplotlib.pyplot as plt


def kotte_nlae(y, p):
    K1ac, K3fdp, L3fdp, K3pep, K2pep, vemax, Kefdp, ne, d, V4max, k1cat, V3max, V2max, ac = p

    flux_1 = k1cat * y[2] * ac / (ac + K1ac)
    flux_2 = vemax * (1 - 1 / (1 + (Kefdp / y[1]) ** ne))
    # convenience kinetics for flux 3
    fdp_sat = y[1] / K3fdp
    pep_sat = y[0] / K3pep
    nr_3 = V3max * fdp_sat
    dr_3 = 1 + fdp_sat
    regulation_activate = 1 / (1 + 1 / pep_sat)
    # regulation_inhibition = 1/(1 + pep_sat) # for future reference
    flux_3 = regulation_activate * nr_3 / dr_3
    flux_4 = V2max * y[0] / (y[0] + K2pep)
    flux_5 = V4max * y[0]
    flux_6 = d * y[2]
    flux = np.hstack((flux_1, flux_2, flux_3, flux_4, flux_5, flux_6))
    # flux = kotte_ck_flux(y, p)

    yd_pep = flux[0] - flux[3] - flux[4]
    yd_fdp = flux[3] - flux[2]
    yd_e = flux[1] - flux[5]

    return np.hstack((yd_pep, yd_fdp, yd_e))


def solve_nlae(x0, model_parameters):
    x_sol, info = opt.fsolve(kotte_nlae, x0, model_parameters, xtol=1e-6)
    return x_sol


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
    x = solve_nlae(y0, parameter_values)
    f, ax = plt.subplots(1, 1)
    # ax.plot()