import numpy as np
import os.path
import csv
from identifiability_analysis import truncate_values
from identifiability_analysis import multi_sample_ident_fun
from identifiability_analysis import get_ident_value

# K1ac, K3fdp, L3fdp, K3pep, K2pep, vemax, Kefdp, ne, d, V4max, k1cat, V3max, V2max, ac
# def_par_val = np.array([.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, .1])
# def_par_val = {"K1ac": np.array([.1]), "K3fdp": np.array([.1]), "L3fdp": np.array([4e6]), "K3pep": np.array([.1]),
#                "K2pep": np.array([.3]),
#                "vemax": np.array([1.1]), "Kefdp": np.array([.45]), "ne": np.array([2]), "d": np.array([.25]),
#                "V4max": np.array([.2]),
#                "k1cat": np.array([1]), "V3max": np.array([1]), "V2max": np.array([1]),
#                "ac": np.array([.1])}


def kotte_true_parameter_values(flux_based=0, flux_name=(), flux_choice_id=0, parameter_id=()):
    if flux_based:
        alter_parameter_value = {"flux1": [{"K1ac (enz)": np.array([.1]), "k1cat": np.array([1]),
                                           "V1max": np.array([1]), "K1ac (no enz)": np.array([.1])},
                                           {"V1max": np.array([1]), "K1ac (no enz)": np.array([.1])},
                                           {"K1ac (enz)": np.array([.1]), "k1cat": np.array([1])}],
                                 "flux2": [{"K2pep": np.array([.3]), "V2max": np.array([1])}],
                                 "flux3": [{"K3fdp (1)": np.array([.1]), "K3pep (1)": np.array([.1]),
                                            "V3max (1)": np.array([1]), "K3fdp (2)": np.array([.1]),
                                            "K3pep (2)": np.array([.1]), "V3max (2)": np.array([1])},
                                           {"K3fdp (1)": np.array([.1]), "K3pep (1)": np.array([.1]),
                                            "V3max (1)": np.array([1])},
                                           {"K3fdp (2)": np.array([.1]), "K3pep (2)": np.array([.1]),
                                            "V3max (2)": np.array([1])}],
                                 "flux5": [{"vemax (1)": np.array([1.1]), "Kefdp (1)": np.array([.45]),
                                            "Kefdp (2)": np.array([.45])},
                                           {"vemax (1)": np.array([1.1]), "Kefdp (1)": np.array([.45])},
                                           {"vemax (1)": np.array([1.1]), "Kefdp (2)": np.array([.45])}]}
        try:
            parameter_value = [alter_parameter_value[name][choice_id][id]
                               for name, choice_id, id in zip(flux_name, flux_choice_id, parameter_id)]
        except TypeError:
            parameter_value = alter_parameter_value[flux_name][flux_choice_id][parameter_id]
        return parameter_value
    else:
        default_parameter_value = {"K1ac": np.array([.1]), "K3fdp": np.array([.1]), "L3fdp": np.array([4e6]),
                                   "K3pep": np.array([.1]), "K2pep": np.array([.3]), "vemax": np.array([1.1]),
                                   "Kefdp": np.array([.45]), "ne": np.array([2]), "d": np.array([.25]),
                                   "V4max": np.array([.2]), "k1cat": np.array([1]), "V3max": np.array([1]),
                                   "V2max": np.array([1]), "ac": np.array([.1])}
        return default_parameter_value


def kotte_ck_flux(y, p={}):
    """calculate flux using convenience kinetics"""

    # K1ac, K3fdp, L3fdp, K3pep, K2pep, vemax, Kefdp, ne, d, V4max, k1cat, V3max, V2max, ac = p
    if not p:
        p = kotte_true_parameter_values()
    K1ac = p["K1ac"]
    K3fdp = p["K3fdp"]
    K3pep = p["K3pep"]
    K2pep = p["K2pep"]
    vemax = p["vemax"]
    Kefdp = p["Kefdp"]
    ne = p["ne"]
    d = p["d"]
    V4max = p["V4max"]
    k1cat = p["k1cat"]
    V3max = p["V3max"]
    V2max = p["V2max"]
    ac = p["ac"]

    flux_1 = k1cat * y[2] * ac / (ac + K1ac)
    flux_2 = vemax * (1 - 1 / (1 + (Kefdp / y[1]) ** ne))
    # convenience kinetics for flux 3
    fdp_sat = y[1] / K3fdp
    pep_sat = y[0] / K3pep
    nr_3 = V3max*fdp_sat
    dr_3 = 1 + fdp_sat
    regulation_activate = 1/(1 + 1/pep_sat)
    # regulation_inhibition = 1/(1 + pep_sat) # for future reference
    flux_3 = regulation_activate * nr_3/dr_3
    flux_4 = V2max * y[0] / (y[0] + K2pep)
    flux_5 = V4max * y[0]
    flux_6 = d * y[2]
    all_flux = np.hstack((flux_1, flux_2, flux_3, flux_4, flux_5, flux_6))

    return all_flux


def kotte_ck_ode(t, y, par_val):
    """ode calculation using convenience kinetics for flux 3"""

    flux = kotte_ck_flux(y, par_val)
    yd_pep = flux[0] - flux[3] - flux[4]
    yd_fdp = flux[3] - flux[2]
    yd_e = flux[1] - flux[5]

    return np.hstack((yd_pep, yd_fdp, yd_e))


def kotte_flux(y, p={}):
    """function doc_string"""

    # K1ac, K3fdp, L3fdp, K3pep, K2pep, vemax, Kefdp, ne, d, V4max, k1cat, V3max, V2max, ac = p
    if not p:
        p = kotte_true_parameter_values()
    K1ac = p["K1ac"]
    K3fdp = p["K3fdp"]
    K3pep = p["K3pep"]
    L3fdp = p["L3fdp"]
    K2pep = p["K2pep"]
    vemax = p["vemax"]
    Kefdp = p["Kefdp"]
    ne = p["ne"]
    d = p["d"]
    V4max = p["V4max"]
    k1cat = p["k1cat"]
    V3max = p["V3max"]
    V2max = p["V2max"]
    ac = p["ac"]

    flux_1 = k1cat * y[2]*ac / (ac + K1ac)
    flux_2 = vemax * (1 - 1 / (1 + (Kefdp / y[1])**ne))
    fdp_sat = 1 + y[1]/K3fdp
    pep_sat = 1 + y[0]/K3pep
    flux_3 = V3max * (fdp_sat - 1) * (fdp_sat**3) / (fdp_sat**4 + L3fdp * (pep_sat**(-4)))
    flux_4 = V2max * y[0] / (y[0]+K2pep)
    flux_5 = V4max * y[0]
    flux_6 = d * y[2]
    all_flux = np.hstack((flux_1, flux_2, flux_3, flux_4, flux_5, flux_6))

    return all_flux


def kotte_ode(t, y, par_val):

    # K1ac, K3fdp, L3fdp, K3pep, K2pep, vemax, Kefdp, ne, d, V4max, k1cat, V3max, V2max, ac = \
    #     [.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, .1]
    # par_val = np.vstack((K1ac, K3fdp, L3fdp, K3pep, K2pep, vemax, Kefdp, ne, d, V4max, k1cat, V3max, V2max, ac))

    flux = kotte_flux(y,par_val)
    yd_pep = flux[0] - flux[3] - flux[4]
    yd_fdp = flux[3] - flux[2]
    yd_e = flux[1] - flux[5]

    return np.hstack((yd_pep, yd_fdp, yd_e))


def flux_1_Vmax_ident(experimental_data):
    """v1 identifiability without enzyme"""
    ac1, _, _, _, v11, _, _, _, _, _, \
    ac2, _, _, _, v12, _, _, _, _, _ = list(experimental_data)

    # flux numerator and denominator w/o sympy
    # symbolic expression for flux v1 w/o enzyme concentration data
    v1max_no_enzyme_numerator_value = ac1 * v11 * v12 - ac2 * v11 * v12
    v1max_no_enzyme_denominator_value = -(ac2 * v11 - ac1 * v12)
    k1ac_no_enzyme_numerator_value = ac1 * (ac2 * v11 - ac2 * v12)
    k1ac_no_enzyme_denominator_value = -ac2 * v11 + ac1 * v12

    v1max_no_enzyme_value = v1max_no_enzyme_numerator_value / v1max_no_enzyme_denominator_value
    k1ac_no_enzyme_value = k1ac_no_enzyme_numerator_value / k1ac_no_enzyme_denominator_value
    return [v1max_no_enzyme_numerator_value, v1max_no_enzyme_denominator_value, v1max_no_enzyme_value], \
           [k1ac_no_enzyme_numerator_value, k1ac_no_enzyme_denominator_value, k1ac_no_enzyme_value]


def flux_1_kcat_ident(experimental_data):
    """v1 identifiability with enzyme concentration"""
    ac1, _, _, x31, v11, _, _, _, _, _, \
    ac2, _, _, x32, v12, _, _, _, _, _ = list(experimental_data)

    k1cat_enzyme_numerator_value = - ac1 * v11 * v12 + ac2 * v11 * v12
    k1cat_enzyme_denominator_value = -(ac1 * v12 * x31 - ac2 * v11 * x32)
    k1cat_enzyme_value = k1cat_enzyme_numerator_value / k1cat_enzyme_denominator_value
    k1ac_enzyme_numerator_value = ac1 * (-ac2 * v12 * x31 + ac2 * v11 * x32)
    k1ac_enzyme_denominator_value = ac1 * v12 * x31 - ac2 * v11 * x32
    k1ac_enzyme_value = k1ac_enzyme_numerator_value / k1ac_enzyme_denominator_value
    return [k1cat_enzyme_numerator_value, k1cat_enzyme_denominator_value, k1cat_enzyme_value], \
           [k1ac_enzyme_numerator_value, k1ac_enzyme_denominator_value, k1ac_enzyme_value]


def flux_1_ident_expression(experimental_data):
    """symbolic and lambdify expression for flux 1 denominator from mathematica"""
    # get variable values (w/o sympy directly from experimental data)
    ac1, x11, x21, x31, v11, v21, v31, v41, v51, v61, \
    ac2, x12, x22, x32, v12, v22, v32, v42, v52, v62 = list(experimental_data)

    # flux numerator and denominator w/o sympy
    # symbolic expression for flux v1 w/o enzyme concentration data
    v1max_no_enzyme_numerator_value = ac1 * v11 * v12 - ac2 * v11 * v12
    v1max_no_enzyme_denominator_value = -(ac2 * v11 - ac1 * v12)
    k1ac_no_enzyme_numerator_value = ac1 * (ac2 * v11 - ac2 * v12)
    k1ac_no_enzyme_denominator_value = -ac2 * v11 + ac1 * v12

    v1max_no_enzyme_value = v1max_no_enzyme_numerator_value/v1max_no_enzyme_denominator_value
    k1ac_no_enzyme_value = k1ac_no_enzyme_numerator_value/k1ac_no_enzyme_denominator_value

    # symbolic expression for flux v1 w/ enzyme concentration data
    k1cat_enzyme_numerator_value = - ac1 * v11 * v12 + ac2 * v11 * v12
    k1cat_enzyme_denominator_value = -(ac1 * v12 * x31 - ac2 * v11 * x32)
    k1cat_enzyme_value = k1cat_enzyme_numerator_value/k1cat_enzyme_denominator_value
    k1ac_enzyme_numerator_value = ac1 * (-ac2 * v12 * x31 + ac2 * v11 * x32)
    k1ac_enzyme_denominator_value = ac1 * v12 * x31 - ac2 * v11 * x32
    k1ac_enzyme_value = k1ac_enzyme_numerator_value/k1ac_enzyme_denominator_value

    return [v1max_no_enzyme_numerator_value, v1max_no_enzyme_denominator_value, v1max_no_enzyme_value], \
           [k1ac_no_enzyme_numerator_value, k1ac_no_enzyme_denominator_value, k1ac_no_enzyme_value], \
           [k1cat_enzyme_numerator_value, k1cat_enzyme_denominator_value, k1cat_enzyme_value], \
           [k1ac_enzyme_numerator_value, k1ac_enzyme_denominator_value, k1ac_enzyme_value]


def flux_2_ident_expression(experimental_data):
    """symbolic and lambdify expression for flux 2 denominator from mathematica"""
    # get variable values (w/o sympy directly from experimental data)
    _, x21, _, _, _, v21, _, _, _, _, \
    _, x22, _, _, _, v22, _, _, _, _ = list(experimental_data)

    # symbolic expression for v2
    # V2max
    v2max_numerator_value = -v21 * v22 * x21 + v21 * v22 * x22
    v2max_denominator_value = -(v22 * x21 - v21 * x22)
    # K2pep
    k2pep_numerator_value = x21 * (v21 * x22 - v22 * x22)
    k2pep_denominator_value = v22 * x21 - v21 * x22
    v2max_value = v2max_numerator_value/v2max_denominator_value
    k2pep_value = k2pep_numerator_value/k2pep_denominator_value

    return [v2max_numerator_value, v2max_denominator_value, v2max_value], \
           [k2pep_numerator_value, k2pep_denominator_value, k2pep_value]


def v3_Vmax_value1(experimental_data):
    """V3max (root 1) identifiability expression for v3"""
    _, x11, x21, _, _, _, v31, _, _, _, \
    _, x12, x22, _, _, _, v32, _, _, _, \
    _, x13, x23, _, _, _, v33, _, _, _ = list(experimental_data)

    # V3max
    sqrt_v3max_nr_1 = ((v31 * v33 * x11 * x12 * x21 * x22 - v32 * v33 * x11 * x12 * x21 * x22 -
                        v32 * v33 * x11 * x13 * x21 * x22 + v31 * v33 * x12 * x13 * x21 * x22 +
                        v32 * v33 * x11 * x12 * x21 * x23 - v31 * v32 * x11 * x13 * x21 * x23 +
                        v32 * v33 * x11 * x13 * x21 * x23 - v31 * v32 * x12 * x13 * x21 * x23 -
                        v31 * v33 * x11 * x12 * x22 * x23 + v31 * v32 * x11 * x13 * x22 * x23 +
                        v31 * v32 * x12 * x13 * x22 * x23 - v31 * v33 * x12 * x13 * x22 * x23) ** 2 -
                       4 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                            v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23) *
                       (v31 * v33 * x11 * x12 * x13 * x21 * x22 - v32 * v33 * x11 * x12 * x13 * x21 * x22 -
                        v31 * v32 * x11 * x12 * x13 * x21 * x23 + v32 * v33 * x11 * x12 * x13 * x21 * x23 +
                        v31 * v32 * x11 * x12 * x13 * x22 * x23 - v31 * v33 * x11 * x12 * x13 * x22 * x23))
    sqrt_v3max_nr_1 = truncate_values(sqrt_v3max_nr_1, 10)
    v3max_nr_1_value = - v31 * v32 * v33 * x11 * x12 * x21 + v31 * v32 * v33 * x11 * x13 * x21 + \
                       v31 * v32 * v33 * x11 * x12 * x22 - v31 * v32 * v33 * x12 * x13 * x22 - \
                       v31 * v32 * v33 * x11 * x13 * x23 + v31 * v32 * v33 * x12 * x13 * x23 - \
                       (v31 * v32 * v33 * x12 * x21 * (-v31 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x13 * x21 * x22 -
                                                       v31 * v33 * x12 * x13 * x21 * x22 -
                                                       v32 * v33 * x11 * x12 * x21 * x23 +
                                                       v31 * v32 * x11 * x13 * x21 * x23 -
                                                       v32 * v33 * x11 * x13 * x21 * x23 +
                                                       v31 * v32 * x12 * x13 * x21 * x23 +
                                                       v31 * v33 * x11 * x12 * x22 * x23 -
                                                       v31 * v32 * x11 * x13 * x22 * x23 -
                                                       v31 * v32 * x12 * x13 * x22 * x23 +
                                                       v31 * v33 * x12 * x13 * x22 * x23 -
                                                       np.sqrt(sqrt_v3max_nr_1))) / \
                       (2 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                             v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23)) + \
                       (v31 * v32 * v33 * x13 * x21 * (-v31 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x13 * x21 * x22 -
                                                       v31 * v33 * x12 * x13 * x21 * x22 -
                                                       v32 * v33 * x11 * x12 * x21 * x23 +
                                                       v31 * v32 * x11 * x13 * x21 * x23 -
                                                       v32 * v33 * x11 * x13 * x21 * x23 +
                                                       v31 * v32 * x12 * x13 * x21 * x23 +
                                                       v31 * v33 * x11 * x12 * x22 * x23 -
                                                       v31 * v32 * x11 * x13 * x22 * x23 -
                                                       v31 * v32 * x12 * x13 * x22 * x23 +
                                                       v31 * v33 * x12 * x13 * x22 * x23 - np.sqrt(sqrt_v3max_nr_1))) / \
                       (2 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                             v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23)) + \
                       (v31 * v32 * v33 * x11 * x22 * (-v31 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x13 * x21 * x22 -
                                                       v31 * v33 * x12 * x13 * x21 * x22 -
                                                       v32 * v33 * x11 * x12 * x21 * x23 +
                                                       v31 * v32 * x11 * x13 * x21 * x23 -
                                                       v32 * v33 * x11 * x13 * x21 * x23 +
                                                       v31 * v32 * x12 * x13 * x21 * x23 +
                                                       v31 * v33 * x11 * x12 * x22 * x23 -
                                                       v31 * v32 * x11 * x13 * x22 * x23 -
                                                       v31 * v32 * x12 * x13 * x22 * x23 +
                                                       v31 * v33 * x12 * x13 * x22 * x23 - np.sqrt(sqrt_v3max_nr_1))) / \
                       (2 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                             v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23)) - \
                       (v31 * v32 * v33 * x13 * x22 * (-v31 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x13 * x21 * x22 -
                                                       v31 * v33 * x12 * x13 * x21 * x22 -
                                                       v32 * v33 * x11 * x12 * x21 * x23 +
                                                       v31 * v32 * x11 * x13 * x21 * x23 -
                                                       v32 * v33 * x11 * x13 * x21 * x23 +
                                                       v31 * v32 * x12 * x13 * x21 * x23 +
                                                       v31 * v33 * x11 * x12 * x22 * x23 -
                                                       v31 * v32 * x11 * x13 * x22 * x23 -
                                                       v31 * v32 * x12 * x13 * x22 * x23 +
                                                       v31 * v33 * x12 * x13 * x22 * x23 - np.sqrt(sqrt_v3max_nr_1))) / \
                       (2 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                             v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23)) - \
                       (v31 * v32 * v33 * x11 * x23 * (-v31 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x13 * x21 * x22 -
                                                       v31 * v33 * x12 * x13 * x21 * x22 -
                                                       v32 * v33 * x11 * x12 * x21 * x23 +
                                                       v31 * v32 * x11 * x13 * x21 * x23 -
                                                       v32 * v33 * x11 * x13 * x21 * x23 +
                                                       v31 * v32 * x12 * x13 * x21 * x23 +
                                                       v31 * v33 * x11 * x12 * x22 * x23 -
                                                       v31 * v32 * x11 * x13 * x22 * x23 -
                                                       v31 * v32 * x12 * x13 * x22 * x23 +
                                                       v31 * v33 * x12 * x13 * x22 * x23 - np.sqrt(sqrt_v3max_nr_1))) / \
                       (2 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                             v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23)) + \
                       (v31 * v32 * v33 * x12 * x23 * (-v31 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x13 * x21 * x22 -
                                                       v31 * v33 * x12 * x13 * x21 * x22 -
                                                       v32 * v33 * x11 * x12 * x21 * x23 +
                                                       v31 * v32 * x11 * x13 * x21 * x23 -
                                                       v32 * v33 * x11 * x13 * x21 * x23 +
                                                       v31 * v32 * x12 * x13 * x21 * x23 +
                                                       v31 * v33 * x11 * x12 * x22 * x23 -
                                                       v31 * v32 * x11 * x13 * x22 * x23 -
                                                       v31 * v32 * x12 * x13 * x22 * x23 +
                                                       v31 * v33 * x12 * x13 * x22 * x23 - np.sqrt(sqrt_v3max_nr_1))) / \
                       (2 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                             v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23))
    v3max_dr_1_value = -v32 * v33 * x11 * x12 * x21 + v32 * v33 * x11 * x13 * x21 + v31 * v33 * x11 * x12 * x22 - \
                       v31 * v33 * x12 * x13 * x22 - v31 * v32 * x11 * x13 * x23 + v31 * v32 * x12 * x13 * x23
    v3max_1_value = v3max_nr_1_value / v3max_dr_1_value
    return [v3max_nr_1_value, v3max_dr_1_value, v3max_1_value]


def v3_Vmax_value2(experimental_data):
    """V3max (root 2) identifiability expression for v3"""
    _, x11, x21, _, _, _, v31, _, _, _, \
    _, x12, x22, _, _, _, v32, _, _, _, \
    _, x13, x23, _, _, _, v33, _, _, _ = list(experimental_data)

    # v3max = second solution
    sqrt_v3max_nr_2 = ((v31 * v33 * x11 * x12 * x21 * x22 - v32 * v33 * x11 * x12 * x21 * x22 -
                        v32 * v33 * x11 * x13 * x21 * x22 + v31 * v33 * x12 * x13 * x21 * x22 +
                        v32 * v33 * x11 * x12 * x21 * x23 - v31 * v32 * x11 * x13 * x21 * x23 +
                        v32 * v33 * x11 * x13 * x21 * x23 - v31 * v32 * x12 * x13 * x21 * x23 -
                        v31 * v33 * x11 * x12 * x22 * x23 + v31 * v32 * x11 * x13 * x22 * x23 +
                        v31 * v32 * x12 * x13 * x22 * x23 - v31 * v33 * x12 * x13 * x22 * x23) ** 2 -
                       4 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                            v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23) *
                       (v31 * v33 * x11 * x12 * x13 * x21 * x22 - v32 * v33 * x11 * x12 * x13 * x21 * x22 -
                        v31 * v32 * x11 * x12 * x13 * x21 * x23 + v32 * v33 * x11 * x12 * x13 * x21 * x23 +
                        v31 * v32 * x11 * x12 * x13 * x22 * x23 - v31 * v33 * x11 * x12 * x13 * x22 * x23))
    sqrt_v3max_nr_2 = truncate_values(sqrt_v3max_nr_2, 10)
    v3max_nr_2_value = -v31 * v32 * v33 * x11 * x12 * x21 + v31 * v32 * v33 * x11 * x13 * x21 + \
                       v31 * v32 * v33 * x11 * x12 * x22 - v31 * v32 * v33 * x12 * x13 * x22 - \
                       v31 * v32 * v33 * x11 * x13 * x23 + v31 * v32 * v33 * x12 * x13 * x23 - \
                       (v31 * v32 * v33 * x12 * x21 * (-v31 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x13 * x21 * x22 -
                                                       v31 * v33 * x12 * x13 * x21 * x22 -
                                                       v32 * v33 * x11 * x12 * x21 * x23 +
                                                       v31 * v32 * x11 * x13 * x21 * x23 -
                                                       v32 * v33 * x11 * x13 * x21 * x23 +
                                                       v31 * v32 * x12 * x13 * x21 * x23 +
                                                       v31 * v33 * x11 * x12 * x22 * x23 -
                                                       v31 * v32 * x11 * x13 * x22 * x23 -
                                                       v31 * v32 * x12 * x13 * x22 * x23 +
                                                       v31 * v33 * x12 * x13 * x22 * x23 + np.sqrt(sqrt_v3max_nr_2))) / \
                       (2 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                             v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23)) + \
                       (v31 * v32 * v33 * x13 * x21 * (-v31 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x13 * x21 * x22 -
                                                       v31 * v33 * x12 * x13 * x21 * x22 -
                                                       v32 * v33 * x11 * x12 * x21 * x23 +
                                                       v31 * v32 * x11 * x13 * x21 * x23 -
                                                       v32 * v33 * x11 * x13 * x21 * x23 +
                                                       v31 * v32 * x12 * x13 * x21 * x23 +
                                                       v31 * v33 * x11 * x12 * x22 * x23 -
                                                       v31 * v32 * x11 * x13 * x22 * x23 -
                                                       v31 * v32 * x12 * x13 * x22 * x23 +
                                                       v31 * v33 * x12 * x13 * x22 * x23 + np.sqrt(sqrt_v3max_nr_2))) / \
                       (2 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                             v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23)) + \
                       (v31 * v32 * v33 * x11 * x22 * (-v31 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x13 * x21 * x22 -
                                                       v31 * v33 * x12 * x13 * x21 * x22 -
                                                       v32 * v33 * x11 * x12 * x21 * x23 +
                                                       v31 * v32 * x11 * x13 * x21 * x23 -
                                                       v32 * v33 * x11 * x13 * x21 * x23 +
                                                       v31 * v32 * x12 * x13 * x21 * x23 +
                                                       v31 * v33 * x11 * x12 * x22 * x23 -
                                                       v31 * v32 * x11 * x13 * x22 * x23 -
                                                       v31 * v32 * x12 * x13 * x22 * x23 +
                                                       v31 * v33 * x12 * x13 * x22 * x23 + np.sqrt(sqrt_v3max_nr_2))) / \
                       (2 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                             v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23)) - \
                       (v31 * v32 * v33 * x13 * x22 * (-v31 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x13 * x21 * x22 -
                                                       v31 * v33 * x12 * x13 * x21 * x22 -
                                                       v32 * v33 * x11 * x12 * x21 * x23 +
                                                       v31 * v32 * x11 * x13 * x21 * x23 -
                                                       v32 * v33 * x11 * x13 * x21 * x23 +
                                                       v31 * v32 * x12 * x13 * x21 * x23 +
                                                       v31 * v33 * x11 * x12 * x22 * x23 -
                                                       v31 * v32 * x11 * x13 * x22 * x23 -
                                                       v31 * v32 * x12 * x13 * x22 * x23 +
                                                       v31 * v33 * x12 * x13 * x22 * x23 + np.sqrt(sqrt_v3max_nr_2))) / \
                       (2 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                             v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23)) - \
                       (v31 * v32 * v33 * x11 * x23 * (-v31 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x13 * x21 * x22 -
                                                       v31 * v33 * x12 * x13 * x21 * x22 -
                                                       v32 * v33 * x11 * x12 * x21 * x23 +
                                                       v31 * v32 * x11 * x13 * x21 * x23 -
                                                       v32 * v33 * x11 * x13 * x21 * x23 +
                                                       v31 * v32 * x12 * x13 * x21 * x23 +
                                                       v31 * v33 * x11 * x12 * x22 * x23 -
                                                       v31 * v32 * x11 * x13 * x22 * x23 -
                                                       v31 * v32 * x12 * x13 * x22 * x23 +
                                                       v31 * v33 * x12 * x13 * x22 * x23 + np.sqrt(sqrt_v3max_nr_2))) / \
                       (2 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                             v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23)) + \
                       (v31 * v32 * v33 * x12 * x23 * (-v31 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x13 * x21 * x22 -
                                                       v31 * v33 * x12 * x13 * x21 * x22 -
                                                       v32 * v33 * x11 * x12 * x21 * x23 +
                                                       v31 * v32 * x11 * x13 * x21 * x23 -
                                                       v32 * v33 * x11 * x13 * x21 * x23 +
                                                       v31 * v32 * x12 * x13 * x21 * x23 +
                                                       v31 * v33 * x11 * x12 * x22 * x23 -
                                                       v31 * v32 * x11 * x13 * x22 * x23 -
                                                       v31 * v32 * x12 * x13 * x22 * x23 +
                                                       v31 * v33 * x12 * x13 * x22 * x23 + np.sqrt(sqrt_v3max_nr_2))) / \
                       (2 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                             v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23))
    v3max_dr_2_value = -v32 * v33 * x11 * x12 * x21 + v32 * v33 * x11 * x13 * x21 + v31 * v33 * x11 * x12 * x22 - \
                       v31 * v33 * x12 * x13 * x22 - v31 * v32 * x11 * x13 * x23 + v31 * v32 * x12 * x13 * x23
    v3max_2_value = v3max_nr_2_value / v3max_dr_2_value

    return [v3max_nr_2_value, v3max_dr_2_value, v3max_2_value]


def v3_K3fdp_value1(experimental_data):
    """K3fdp (root 1) identifiability expression for v3"""
    _, x11, x21, _, _, _, v31, _, _, _, \
    _, x12, x22, _, _, _, v32, _, _, _, \
    _, x13, x23, _, _, _, v33, _, _, _ = list(experimental_data)

    # K3fdp
    sq_nr_1_k3fdp = ((v31 * v33 * x11 * x12 * x21 * x22 - v32 * v33 * x11 * x12 * x21 * x22 -
                      v32 * v33 * x11 * x13 * x21 * x22 + v31 * v33 * x12 * x13 * x21 * x22 +
                      v32 * v33 * x11 * x12 * x21 * x23 - v31 * v32 * x11 * x13 * x21 * x23 +
                      v32 * v33 * x11 * x13 * x21 * x23 - v31 * v32 * x12 * x13 * x21 * x23 -
                      v31 * v33 * x11 * x12 * x22 * x23 + v31 * v32 * x11 * x13 * x22 * x23 +
                      v31 * v32 * x12 * x13 * x22 * x23 - v31 * v33 * x12 * x13 * x22 * x23) ** 2 -
                     4 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                          v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23) *
                     (v31 * v33 * x11 * x12 * x13 * x21 * x22 - v32 * v33 * x11 * x12 * x13 * x21 * x22 -
                      v31 * v32 * x11 * x12 * x13 * x21 * x23 + v32 * v33 * x11 * x12 * x13 * x21 * x23 +
                      v31 * v32 * x11 * x12 * x13 * x22 * x23 - v31 * v33 * x11 * x12 * x13 * x22 * x23))
    sq_nr_1_k3fdp = truncate_values(sq_nr_1_k3fdp, 10)
    k3fdp_nr_1_value = -v31 * v33 * x11 * x12 * x21 * x22 + v32 * v33 * x11 * x12 * x21 * x22 + \
                       v31 * v32 * x11 * x13 * x21 * x23 - v32 * v33 * x11 * x13 * x21 * x23 - \
                       v31 * v32 * x12 * x13 * x22 * x23 + v31 * v33 * x12 * x13 * x22 * x23 + \
                       (v32 * v33 * x11 * x21 * x22 * (-v31 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x13 * x21 * x22 -
                                                       v31 * v33 * x12 * x13 * x21 * x22 -
                                                       v32 * v33 * x11 * x12 * x21 * x23 +
                                                       v31 * v32 * x11 * x13 * x21 * x23 -
                                                       v32 * v33 * x11 * x13 * x21 * x23 +
                                                       v31 * v32 * x12 * x13 * x21 * x23 +
                                                       v31 * v33 * x11 * x12 * x22 * x23 -
                                                       v31 * v32 * x11 * x13 * x22 * x23 -
                                                       v31 * v32 * x12 * x13 * x22 * x23 +
                                                       v31 * v33 * x12 * x13 * x22 * x23 - np.sqrt(sq_nr_1_k3fdp))) / \
                       (2 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                             v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23)) - \
                       (v31 * v33 * x12 * x21 * x22 * (-v31 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x13 * x21 * x22 -
                                                       v31 * v33 * x12 * x13 * x21 * x22 -
                                                       v32 * v33 * x11 * x12 * x21 * x23 +
                                                       v31 * v32 * x11 * x13 * x21 * x23 -
                                                       v32 * v33 * x11 * x13 * x21 * x23 +
                                                       v31 * v32 * x12 * x13 * x21 * x23 +
                                                       v31 * v33 * x11 * x12 * x22 * x23 -
                                                       v31 * v32 * x11 * x13 * x22 * x23 -
                                                       v31 * v32 * x12 * x13 * x22 * x23 +
                                                       v31 * v33 * x12 * x13 * x22 * x23 - np.sqrt(sq_nr_1_k3fdp))) / \
                       (2 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                             v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23)) - \
                       (v32 * v33 * x11 * x21 * x23 * (-v31 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x13 * x21 * x22 -
                                                       v31 * v33 * x12 * x13 * x21 * x22 -
                                                       v32 * v33 * x11 * x12 * x21 * x23 +
                                                       v31 * v32 * x11 * x13 * x21 * x23 -
                                                       v32 * v33 * x11 * x13 * x21 * x23 +
                                                       v31 * v32 * x12 * x13 * x21 * x23 +
                                                       v31 * v33 * x11 * x12 * x22 * x23 -
                                                       v31 * v32 * x11 * x13 * x22 * x23 -
                                                       v31 * v32 * x12 * x13 * x22 * x23 +
                                                       v31 * v33 * x12 * x13 * x22 * x23 - np.sqrt(sq_nr_1_k3fdp))) / \
                       (2 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                             v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23)) + \
                       (v31 * v32 * x13 * x21 * x23 * (-v31 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x13 * x21 * x22 -
                                                       v31 * v33 * x12 * x13 * x21 * x22 -
                                                       v32 * v33 * x11 * x12 * x21 * x23 +
                                                       v31 * v32 * x11 * x13 * x21 * x23 -
                                                       v32 * v33 * x11 * x13 * x21 * x23 +
                                                       v31 * v32 * x12 * x13 * x21 * x23 +
                                                       v31 * v33 * x11 * x12 * x22 * x23 -
                                                       v31 * v32 * x11 * x13 * x22 * x23 -
                                                       v31 * v32 * x12 * x13 * x22 * x23 +
                                                       v31 * v33 * x12 * x13 * x22 * x23 - np.sqrt(sq_nr_1_k3fdp))) / \
                       (2 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                             v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23)) + \
                       (v31 * v33 * x12 * x22 * x23 * (-v31 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x13 * x21 * x22 -
                                                       v31 * v33 * x12 * x13 * x21 * x22 -
                                                       v32 * v33 * x11 * x12 * x21 * x23 +
                                                       v31 * v32 * x11 * x13 * x21 * x23 -
                                                       v32 * v33 * x11 * x13 * x21 * x23 +
                                                       v31 * v32 * x12 * x13 * x21 * x23 +
                                                       v31 * v33 * x11 * x12 * x22 * x23 -
                                                       v31 * v32 * x11 * x13 * x22 * x23 -
                                                       v31 * v32 * x12 * x13 * x22 * x23 +
                                                       v31 * v33 * x12 * x13 * x22 * x23 - np.sqrt(sq_nr_1_k3fdp))) / \
                       (2 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                             v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23)) - \
                       (v31 * v32 * x13 * x22 * x23 * (-v31 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x13 * x21 * x22 -
                                                       v31 * v33 * x12 * x13 * x21 * x22 -
                                                       v32 * v33 * x11 * x12 * x21 * x23 +
                                                       v31 * v32 * x11 * x13 * x21 * x23 -
                                                       v32 * v33 * x11 * x13 * x21 * x23 +
                                                       v31 * v32 * x12 * x13 * x21 * x23 +
                                                       v31 * v33 * x11 * x12 * x22 * x23 -
                                                       v31 * v32 * x11 * x13 * x22 * x23 -
                                                       v31 * v32 * x12 * x13 * x22 * x23 +
                                                       v31 * v33 * x12 * x13 * x22 * x23 - np.sqrt(sq_nr_1_k3fdp))) / \
                       (2 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                             v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23))
    k3fdp_dr_1_value = -v32 * v33 * x11 * x12 * x21 + v32 * v33 * x11 * x13 * x21 + v31 * v33 * x11 * x12 * x22 - \
                       v31 * v33 * x12 * x13 * x22 - v31 * v32 * x11 * x13 * x23 + v31 * v32 * x12 * x13 * x23
    k3fdp_1_value = k3fdp_nr_1_value / k3fdp_dr_1_value
    return [k3fdp_nr_1_value, k3fdp_dr_1_value, k3fdp_1_value]


def v3_K3fdp_value2(experimental_data):
    """K3fdp (root 2) identifiability expression for v3"""
    _, x11, x21, _, _, _, v31, _, _, _, \
    _, x12, x22, _, _, _, v32, _, _, _, \
    _, x13, x23, _, _, _, v33, _, _, _ = list(experimental_data)

    # K3fdp 2
    sq_k3fdp_nr_2 = ((v31 * v33 * x11 * x12 * x21 * x22 - v32 * v33 * x11 * x12 * x21 * x22 -
                      v32 * v33 * x11 * x13 * x21 * x22 + v31 * v33 * x12 * x13 * x21 * x22 +
                      v32 * v33 * x11 * x12 * x21 * x23 - v31 * v32 * x11 * x13 * x21 * x23 +
                      v32 * v33 * x11 * x13 * x21 * x23 - v31 * v32 * x12 * x13 * x21 * x23 -
                      v31 * v33 * x11 * x12 * x22 * x23 + v31 * v32 * x11 * x13 * x22 * x23 +
                      v31 * v32 * x12 * x13 * x22 * x23 - v31 * v33 * x12 * x13 * x22 * x23) ** 2 -
                     4 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                          v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23) *
                     (v31 * v33 * x11 * x12 * x13 * x21 * x22 - v32 * v33 * x11 * x12 * x13 * x21 * x22 -
                      v31 * v32 * x11 * x12 * x13 * x21 * x23 + v32 * v33 * x11 * x12 * x13 * x21 * x23 +
                      v31 * v32 * x11 * x12 * x13 * x22 * x23 - v31 * v33 * x11 * x12 * x13 * x22 * x23))
    sq_k3fdp_nr_2 = truncate_values(sq_k3fdp_nr_2, 10)
    k3fdp_nr_2_value = -v31 * v33 * x11 * x12 * x21 * x22 + v32 * v33 * x11 * x12 * x21 * x22 + \
                       v31 * v32 * x11 * x13 * x21 * x23 - v32 * v33 * x11 * x13 * x21 * x23 - \
                       v31 * v32 * x12 * x13 * x22 * x23 + v31 * v33 * x12 * x13 * x22 * x23 + \
                       (v32 * v33 * x11 * x21 * x22 * (-v31 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x13 * x21 * x22 -
                                                       v31 * v33 * x12 * x13 * x21 * x22 -
                                                       v32 * v33 * x11 * x12 * x21 * x23 +
                                                       v31 * v32 * x11 * x13 * x21 * x23 -
                                                       v32 * v33 * x11 * x13 * x21 * x23 +
                                                       v31 * v32 * x12 * x13 * x21 * x23 +
                                                       v31 * v33 * x11 * x12 * x22 * x23 -
                                                       v31 * v32 * x11 * x13 * x22 * x23 -
                                                       v31 * v32 * x12 * x13 * x22 * x23 +
                                                       v31 * v33 * x12 * x13 * x22 * x23 + np.sqrt(sq_k3fdp_nr_2))) / \
                       (2 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                             v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23)) - \
                       (v31 * v33 * x12 * x21 * x22 * (-v31 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x13 * x21 * x22 -
                                                       v31 * v33 * x12 * x13 * x21 * x22 -
                                                       v32 * v33 * x11 * x12 * x21 * x23 +
                                                       v31 * v32 * x11 * x13 * x21 * x23 -
                                                       v32 * v33 * x11 * x13 * x21 * x23 +
                                                       v31 * v32 * x12 * x13 * x21 * x23 +
                                                       v31 * v33 * x11 * x12 * x22 * x23 -
                                                       v31 * v32 * x11 * x13 * x22 * x23 -
                                                       v31 * v32 * x12 * x13 * x22 * x23 +
                                                       v31 * v33 * x12 * x13 * x22 * x23 + np.sqrt(sq_k3fdp_nr_2))) / \
                       (2 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                             v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23)) - \
                       (v32 * v33 * x11 * x21 * x23 * (-v31 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x13 * x21 * x22 -
                                                       v31 * v33 * x12 * x13 * x21 * x22 -
                                                       v32 * v33 * x11 * x12 * x21 * x23 +
                                                       v31 * v32 * x11 * x13 * x21 * x23 -
                                                       v32 * v33 * x11 * x13 * x21 * x23 +
                                                       v31 * v32 * x12 * x13 * x21 * x23 +
                                                       v31 * v33 * x11 * x12 * x22 * x23 -
                                                       v31 * v32 * x11 * x13 * x22 * x23 -
                                                       v31 * v32 * x12 * x13 * x22 * x23 +
                                                       v31 * v33 * x12 * x13 * x22 * x23 + np.sqrt(sq_k3fdp_nr_2))) / \
                       (2 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                             v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23)) + \
                       (v31 * v32 * x13 * x21 * x23 * (-v31 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x13 * x21 * x22 -
                                                       v31 * v33 * x12 * x13 * x21 * x22 -
                                                       v32 * v33 * x11 * x12 * x21 * x23 +
                                                       v31 * v32 * x11 * x13 * x21 * x23 -
                                                       v32 * v33 * x11 * x13 * x21 * x23 +
                                                       v31 * v32 * x12 * x13 * x21 * x23 +
                                                       v31 * v33 * x11 * x12 * x22 * x23 -
                                                       v31 * v32 * x11 * x13 * x22 * x23 -
                                                       v31 * v32 * x12 * x13 * x22 * x23 +
                                                       v31 * v33 * x12 * x13 * x22 * x23 + np.sqrt(sq_k3fdp_nr_2))) / \
                       (2 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                             v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23)) + \
                       (v31 * v33 * x12 * x22 * x23 * (-v31 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x13 * x21 * x22 -
                                                       v31 * v33 * x12 * x13 * x21 * x22 -
                                                       v32 * v33 * x11 * x12 * x21 * x23 +
                                                       v31 * v32 * x11 * x13 * x21 * x23 -
                                                       v32 * v33 * x11 * x13 * x21 * x23 +
                                                       v31 * v32 * x12 * x13 * x21 * x23 +
                                                       v31 * v33 * x11 * x12 * x22 * x23 -
                                                       v31 * v32 * x11 * x13 * x22 * x23 -
                                                       v31 * v32 * x12 * x13 * x22 * x23 +
                                                       v31 * v33 * x12 * x13 * x22 * x23 + np.sqrt(sq_k3fdp_nr_2))) / \
                       (2 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                             v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23)) - \
                       (v31 * v32 * x13 * x22 * x23 * (-v31 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x12 * x21 * x22 +
                                                       v32 * v33 * x11 * x13 * x21 * x22 -
                                                       v31 * v33 * x12 * x13 * x21 * x22 -
                                                       v32 * v33 * x11 * x12 * x21 * x23 +
                                                       v31 * v32 * x11 * x13 * x21 * x23 -
                                                       v32 * v33 * x11 * x13 * x21 * x23 +
                                                       v31 * v32 * x12 * x13 * x21 * x23 +
                                                       v31 * v33 * x11 * x12 * x22 * x23 -
                                                       v31 * v32 * x11 * x13 * x22 * x23 -
                                                       v31 * v32 * x12 * x13 * x22 * x23 +
                                                       v31 * v33 * x12 * x13 * x22 * x23 + np.sqrt(sq_k3fdp_nr_2))) / \
                       (2 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                             v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23))
    k3fdp_dr_2_value = -v32 * v33 * x11 * x12 * x21 + v32 * v33 * x11 * x13 * x21 + v31 * v33 * x11 * x12 * x22 - \
                       v31 * v33 * x12 * x13 * x22 - v31 * v32 * x11 * x13 * x23 + v31 * v32 * x12 * x13 * x23
    k3fdp_2_value = k3fdp_nr_2_value / k3fdp_dr_2_value
    return [k3fdp_nr_2_value, k3fdp_dr_2_value, k3fdp_2_value]


def v3_K3pep_value1(experimental_data):
    """K3pep (root 1) identifiability expression for v3"""
    _, x11, x21, _, _, _, v31, _, _, _, \
    _, x12, x22, _, _, _, v32, _, _, _, \
    _, x13, x23, _, _, _, v33, _, _, _ = list(experimental_data)

    # K3pep
    sq_k3pep_nr_1 = ((v31 * v33 * x11 * x12 * x21 * x22 - v32 * v33 * x11 * x12 * x21 * x22 -
                      v32 * v33 * x11 * x13 * x21 * x22 + v31 * v33 * x12 * x13 * x21 * x22 +
                      v32 * v33 * x11 * x12 * x21 * x23 - v31 * v32 * x11 * x13 * x21 * x23 +
                      v32 * v33 * x11 * x13 * x21 * x23 - v31 * v32 * x12 * x13 * x21 * x23 -
                      v31 * v33 * x11 * x12 * x22 * x23 + v31 * v32 * x11 * x13 * x22 * x23 +
                      v31 * v32 * x12 * x13 * x22 * x23 - v31 * v33 * x12 * x13 * x22 * x23) ** 2 -
                     4 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                          v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23) *
                     (v31 * v33 * x11 * x12 * x13 * x21 * x22 - v32 * v33 * x11 * x12 * x13 * x21 * x22 -
                      v31 * v32 * x11 * x12 * x13 * x21 * x23 + v32 * v33 * x11 * x12 * x13 * x21 * x23 +
                      v31 * v32 * x11 * x12 * x13 * x22 * x23 - v31 * v33 * x11 * x12 * x13 * x22 * x23))
    sq_k3pep_nr_1 = truncate_values(sq_k3pep_nr_1, 10)
    k3pep_nr_1_value = -v31 * v33 * x11 * x12 * x21 * x22 + v32 * v33 * x11 * x12 * x21 * x22 + \
                       v32 * v33 * x11 * x13 * x21 * x22 - v31 * v33 * x12 * x13 * x21 * x22 - \
                       v32 * v33 * x11 * x12 * x21 * x23 + v31 * v32 * x11 * x13 * x21 * x23 - \
                       v32 * v33 * x11 * x13 * x21 * x23 + v31 * v32 * x12 * x13 * x21 * x23 + \
                       v31 * v33 * x11 * x12 * x22 * x23 - v31 * v32 * x11 * x13 * x22 * x23 - \
                       v31 * v32 * x12 * x13 * x22 * x23 + v31 * v33 * x12 * x13 * x22 * x23 - np.sqrt(sq_k3pep_nr_1)
    k3pep_dr_1_value = 2 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                            v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23)
    k3pep_1_value = k3pep_nr_1_value / k3pep_dr_1_value
    return [k3pep_nr_1_value, k3pep_dr_1_value, k3pep_1_value]


def v3_K3pep_value2(experimental_data):
    """K3pep (root 2) identifiability expression for v3"""
    _, x11, x21, _, _, _, v31, _, _, _, \
    _, x12, x22, _, _, _, v32, _, _, _, \
    _, x13, x23, _, _, _, v33, _, _, _ = list(experimental_data)

    # K3pep 2
    sq_k3pep_nr_2 = ((v31 * v33 * x11 * x12 * x21 * x22 - v32 * v33 * x11 * x12 * x21 * x22 -
                      v32 * v33 * x11 * x13 * x21 * x22 + v31 * v33 * x12 * x13 * x21 * x22 +
                      v32 * v33 * x11 * x12 * x21 * x23 - v31 * v32 * x11 * x13 * x21 * x23 +
                      v32 * v33 * x11 * x13 * x21 * x23 - v31 * v32 * x12 * x13 * x21 * x23 -
                      v31 * v33 * x11 * x12 * x22 * x23 + v31 * v32 * x11 * x13 * x22 * x23 +
                      v31 * v32 * x12 * x13 * x22 * x23 - v31 * v33 * x12 * x13 * x22 * x23) ** 2 -
                     4 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                          v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23) *
                     (v31 * v33 * x11 * x12 * x13 * x21 * x22 - v32 * v33 * x11 * x12 * x13 * x21 * x22 -
                      v31 * v32 * x11 * x12 * x13 * x21 * x23 + v32 * v33 * x11 * x12 * x13 * x21 * x23 +
                      v31 * v32 * x11 * x12 * x13 * x22 * x23 - v31 * v33 * x11 * x12 * x13 * x22 * x23))
    sq_k3pep_nr_2 = truncate_values(sq_k3pep_nr_2, 10)
    k3pep_nr_2_value = -v31 * v33 * x11 * x12 * x21 * x22 + v32 * v33 * x11 * x12 * x21 * x22 + \
                       v32 * v33 * x11 * x13 * x21 * x22 - v31 * v33 * x12 * x13 * x21 * x22 - \
                       v32 * v33 * x11 * x12 * x21 * x23 + v31 * v32 * x11 * x13 * x21 * x23 - \
                       v32 * v33 * x11 * x13 * x21 * x23 + v31 * v32 * x12 * x13 * x21 * x23 + \
                       v31 * v33 * x11 * x12 * x22 * x23 - v31 * v32 * x11 * x13 * x22 * x23 - \
                       v31 * v32 * x12 * x13 * x22 * x23 + v31 * v33 * x12 * x13 * x22 * x23 + np.sqrt(sq_k3pep_nr_2)
    k3pep_dr_2_value = 2 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                            v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23)
    k3pep_2_value = k3pep_nr_2_value / k3pep_dr_2_value
    return [k3pep_nr_2_value, k3pep_dr_2_value, k3pep_2_value]


def flux_3_value1_ident(experimental_data):
    """only one set of expressions separated for flux 3"""

    # V3max
    v3max_nr_1_value, v3max_dr_1_value, v3max_1_value = v3_Vmax_value1(experimental_data)

    # K3fdp
    k3fdp_nr_1_value, k3fdp_dr_1_value, k3fdp_1_value = v3_K3fdp_value1(experimental_data)

    # K3pep
    k3pep_nr_1_value, k3pep_dr_1_value, k3pep_1_value = v3_K3pep_value1(experimental_data)

    return [v3max_nr_1_value, v3max_dr_1_value, v3max_1_value], \
           [k3fdp_nr_1_value, k3fdp_dr_1_value, k3fdp_1_value], \
           [k3pep_nr_1_value, k3pep_dr_1_value, k3pep_1_value]


def flux_3_value2_ident(experimental_data):
    """second set of expressions separated for v3"""

    # v3max = second solution
    v3max_nr_2_value, v3max_dr_2_value, v3max_2_value = v3_Vmax_value2(experimental_data)

    # K3fdp 2
    k3fdp_nr_2_value, k3fdp_dr_2_value, k3fdp_2_value = v3_K3fdp_value2(experimental_data)

    # K3pep 2
    k3pep_nr_2_value, k3pep_dr_2_value, k3pep_2_value = v3_K3pep_value2(experimental_data)

    return [v3max_nr_2_value, v3max_dr_2_value, v3max_2_value], \
           [k3fdp_nr_2_value, k3fdp_dr_2_value, k3fdp_2_value], \
           [k3pep_nr_2_value, k3pep_dr_2_value, k3pep_2_value]


def flux_3_ident_expression(experimental_data):
    """symbolic and lambdify expression for flux 3 denominator from mathematica"""

    v3max_value_1, k3fdp_value_1, k3pep_value_1 = flux_3_value1_ident(experimental_data)
    v3max_value_2, k3fdp_value_2, k3pep_value_2 = flux_3_value2_ident(experimental_data)

    return v3max_value_1, k3fdp_value_1, k3pep_value_1, \
           v3max_value_2, k3fdp_value_2, k3pep_value_2


def v3_V3max_var1(experimental_data):
    """v3max values when k3pep is assumed as known"""
    _, x11, x21, _, _, _, v31, _, _, _, \
    _, x12, x22, _, _, _, v32, _, _, _, \
    _, x13, x23, _, _, _, v33, _, _, _ = list(experimental_data)
    K3pep = np.array([.1])
    v3max_nr_1 = -((-(K3pep * v32 + v32 * x12) * (K3pep * v31 * x21 + v31 * x11 * x21) + (K3pep * v31 + v31 * x11) * (
                K3pep * v32 * x22 + v32 * x12 * x22)))
    v3max_dr_1 = (K3pep*v32*x11*x21 + v32*x11*x12*x21 - K3pep*v31*x12*x22 - v31*x11*x12*x22)
    v3max_1 = v3max_nr_1/v3max_dr_1
    return [v3max_nr_1, v3max_dr_1, v3max_1]


def v3_V3max_var2(experimental_data):
    """v3max value when k3fdp value is assumed as known"""
    _, x11, x21, _, _, _, v31, _, _, _, \
    _, x12, x22, _, _, _, v32, _, _, _, \
    _, x13, x23, _, _, _, v33, _, _, _ = list(experimental_data)
    K3fdp = np.array([.1])
    v3max_nr_2 = -((-(K3fdp*v31*x11 + v31*x11*x21)*(K3fdp*v32 + v32*x22) + (K3fdp*v31 + v31*x21)*(K3fdp*v32*x12 + v32*x12*x22)))
    v3max_dr_2 = (K3fdp*v32*x11*x21 - K3fdp*v31*x12*x22 + v32*x11*x21*x22 - v31*x12*x21*x22)
    v3max_2 = v3max_nr_2/v3max_dr_2
    return [v3max_nr_2, v3max_dr_2, v3max_2]


def v3_K3fdp_var1(experimental_data):
    """k3fdp value when k3pep is assumed as known"""
    _, x11, x21, _, _, _, v31, _, _, _, \
    _, x12, x22, _, _, _, v32, _, _, _, \
    _, x13, x23, _, _, _, v33, _, _, _ = list(experimental_data)
    K3pep = np.array([.1])
    k3fdp_nr = (x21*(-K3pep*v32*x11*x22 + K3pep*v31*x12*x22 + v31*x11*x12*x22 - v32*x11*x12*x22))
    k3fdp_dr = (K3pep*v32*x11*x21 + v32*x11*x12*x21 - K3pep*v31*x12*x22 - v31*x11*x12*x22)
    k3fdp = k3fdp_nr/k3fdp_dr
    return [k3fdp_nr, k3fdp_dr, k3fdp]


def v3_K3pep_var2(experimental_data):
    """k3pep value when k3fdp is assumed as known"""
    _, x11, x21, _, _, _, v31, _, _, _, \
    _, x12, x22, _, _, _, v32, _, _, _, \
    _, x13, x23, _, _, _, v33, _, _, _ = list(experimental_data)
    K3fdp = np.array([.1])
    k3pep_nr = (x11*(-K3fdp*v32*x12*x21 + K3fdp*v31*x12*x22 + v31*x12*x21*x22 - v32*x12*x21*x22))
    k3pep_dr = (K3fdp*v32*x11*x21 - K3fdp*v31*x12*x22 + v32*x11*x21*x22 - v31*x12*x21*x22)
    k3pep = k3pep_nr/k3pep_dr
    return [k3pep_nr, k3pep_dr, k3pep]


def flux_3_var1(experimental_data):
    v3max_nr_1, v3max_dr_1, v3max_1 = v3_V3max_var1(experimental_data)
    k3fdp_nr, k3fdp_dr, k3fdp = v3_K3fdp_var1(experimental_data)
    return [v3max_nr_1, v3max_dr_1, v3max_1], [k3fdp_nr, k3fdp_dr, k3fdp]


def flux_3_var2(experimental_data):
    v3max_nr_2, v3max_dr_2, v3max_2 = v3_V3max_var2(experimental_data)
    k3pep_nr, k3pep_dr, k3pep = v3_K3pep_var2(experimental_data)
    return [v3max_nr_2, v3max_dr_2, v3max_2], [k3pep_nr, k3pep_dr, k3pep]


def flux_3_var_1_and_2(experimental_data):
    v3max_nr_1, v3max_dr_1, v3max_1 = v3_V3max_var1(experimental_data)
    k3fdp_nr, k3fdp_dr, k3fdp = v3_K3fdp_var1(experimental_data)
    v3max_nr_2, v3max_dr_2, v3max_2 = v3_V3max_var2(experimental_data)
    k3pep_nr, k3pep_dr, k3pep = v3_K3pep_var2(experimental_data)

    return [v3max_nr_1, v3max_dr_1, v3max_1], [k3fdp_nr, k3fdp_dr, k3fdp], \
           [v3max_nr_2, v3max_dr_2, v3max_2], [k3pep_nr, k3pep_dr, k3pep]


def v5_vemax_value1_ident(experimental_data):
    """both value 1 and value 2 are the same for vemax"""
    _, _, x21, x31, _, _, _, _, v51, v61, \
    _, _, x22, x32, _, _, _, _, v52, v62 = list(experimental_data)

    vemax_nr = v51 * v52 * (x21**2 - x22**2)
    vemax_dr = v51*x21**2 - v52*x22**2
    vemax_value = vemax_nr/vemax_dr

    return [vemax_nr, vemax_dr, vemax_value]


def v5_Kefdp_value1_ident(experimental_data):
    """Kefdp has 2 different values for v5 upon determination of
    closed-form expressions due to the presence of a square root"""
    _, _, x21, x31, _, _, _, _, v51, v61, \
    _, _, x22, x32, _, _, _, _, v52, v62 = list(experimental_data)

    kefdp_nr_value1 = - np.sqrt(-v51*x21**2 + v52*x22**2)
    kefdp_dr_value1 = np.sqrt(v51 - v52)
    kefdp_value1 = kefdp_nr_value1/kefdp_dr_value1

    return [kefdp_nr_value1, kefdp_dr_value1, kefdp_value1]


def v5_Kefdp_value2_ident(experimental_data):
    _, _, x21, x31, _, _, _, _, v51, v61, \
    _, _, x22, x32, _, _, _, _, v52, v62 = list(experimental_data)

    kefdp_nr_value2 = np.sqrt(-v51 * x21 ** 2 + v52 * x22 ** 2)
    kefdp_dr_value2 = np.sqrt(v51 - v52)
    kefdp_value2 = kefdp_nr_value2 / kefdp_dr_value2

    return [kefdp_nr_value2, kefdp_dr_value2, kefdp_value2]


def flux_5_value1_ident(experimental_data):
    """identifiability expression for transcriptional regulatory reaction (v5) with Hill Kinetics
    Kefdp value 1"""
    vemax_nr_1, vemax_dr_1, vemax_1_value = v5_vemax_value1_ident(experimental_data)
    kefdp_nr_1, kefdp_dr_1, kefdp_1_value = v5_Kefdp_value1_ident(experimental_data)

    return [vemax_nr_1, vemax_dr_1, vemax_1_value], [kefdp_nr_1, kefdp_dr_1, kefdp_1_value]


def flux_5_value2_ident(experimental_data):
    """identifiability expression for transcriptional regulatory reaction (v5) with Hill Kinetics
    Kefdp value 2"""
    vemax_nr_1, vemax_dr_1, vemax_1_value = v5_vemax_value1_ident(experimental_data)
    kefdp_nr_2, kefdp_dr_2, kefdp_2_value = v5_Kefdp_value2_ident(experimental_data)

    return [vemax_nr_1, vemax_dr_1, vemax_1_value], [kefdp_nr_2, kefdp_dr_2, kefdp_2_value]


def flux_5_ident_expression(experimental_data):
    """flux 5 (v5) identifiability for both values simultaneously"""
    # value 1
    vemax_nr_1, vemax_dr_1, vemax_1_value = v5_vemax_value1_ident(experimental_data)
    kefdp_nr_1, kefdp_dr_1, kefdp_1_value = v5_Kefdp_value1_ident(experimental_data)

    # value 2
    vemax_nr_2, vemax_dr_2, vemax_2_value = v5_vemax_value1_ident(experimental_data)
    kefdp_nr_2, kefdp_dr_2, kefdp_2_value = v5_Kefdp_value2_ident(experimental_data)

    return [vemax_nr_1, vemax_dr_1, vemax_1_value], [kefdp_nr_1, kefdp_dr_1, kefdp_1_value], \
           [vemax_nr_2, vemax_dr_2, vemax_2_value], [kefdp_nr_2, kefdp_dr_2, kefdp_2_value]


def ident_parameter_name(parameter_id, flux_name=(), flux_choice_id=0):
    parameter_list = ['V1max', 'K1ac (no enz)', 'k1cat', 'K1ac (enz)',
                      'V2max', 'K2pep',
                      'V3max (1)', 'K3fdp (1)', 'K3pep (1)',
                      'V3max (2)', 'K3fdp (2)', 'K3pep (2)',
                      'vemax (1)', 'Kefdp (1)',
                      'vemax(1)', 'Kefdp (2)']
    # flux_parameter_list = {"flux1":['V1max', 'K1ac (no enz)', 'k1cat', 'K1ac (enz)'],
    #                        "flux2":['V2max', 'K2pep'],
    #                        "flux3":['V3max (1)', 'K3fdp (1)', 'K3pep (1)',
    #                                 'V3max (2)', 'K3fdp (2)', 'K3pep (2)']}
    alter_list = {"flux1": [['V1max', 'K1ac (no enz)', 'k1cat', 'K1ac (enz)'],
                            ['V1max', 'K1ac (no enz)'],
                            ['k1cat', 'K1ac (enz)']],
                  "flux2": [['V2max', 'K2pep']],
                  "flux3": [['V3max (1)', 'K3fdp (1)', 'K3pep (1)', 'V3max (2)', 'K3fdp (2)', 'K3pep (2)'],
                            ['V3max (1)', 'K3fdp (1)', 'K3pep (1)'],
                            ['V3max (2)', 'K3fdp (2)', 'K3pep (2)']],
                  "flux4": [],
                  "flux5": [['vemax (1)', 'Kefdp (1)', 'vemax (1)', 'Kefdp (2)'],
                            ['vemax (1)', 'Kefdp (1)'],
                            ['vemax (1)', 'Kefdp (2)']],
                  "flux6": []}
    if flux_name:
        try:
            parameter_name = [alter_list[name][choice_id][id]
                              for name, choice_id, id in zip(flux_name, flux_choice_id, parameter_id)]
        except TypeError:
            parameter_name = alter_list[flux_name][flux_choice_id][parameter_id]
    else:
        try:
            parameter_name = [parameter_list[id] for id in parameter_id]
        except TypeError:
            parameter_name = parameter_list[parameter_id]
    return parameter_name


def flux_based_id(parameter_id):
    parameter_list = ['V1max', 'K1ac (no enz)', 'k1cat', 'K1ac (enz)',
                      'V2max', 'K2pep',
                      'V3max (1)', 'K3fdp (1)', 'K3pep (1)',
                      'V3max (2)', 'K3fdp (2)', 'K3pep (2)',
                      'vemax (1)', 'Kefdp (1)',
                      'vemax (1)', 'Kefdp (2)']
    flux_parameter_list = {"flux1": ['V1max', 'K1ac (no enz)', 'k1cat', 'K1ac (enz)'],
                           "flux2": ['V2max', 'K2pep'],
                           "flux3": ['V3max (1)', 'K3fdp (1)', 'K3pep (1)',
                                     'V3max (2)', 'K3fdp (2)', 'K3pep (2)'],
                           "flux4": [],
                           "flux5": ['vemax (1)', 'Kefdp (1)', 'vemax(1)', 'Kefdp (2)'],
                           "flux6": []}
    try:
        parameter_name = [(id, parameter_list[id]) for id in parameter_id]
    except TypeError:
        parameter_name = (parameter_id, parameter_list[parameter_id])
    chosen_flux_name = []
    chosen_flux_parameter_id = []
    chosen_flux_orig_id = []
    if isinstance(parameter_name, list):
        for orig_id, j_parameter in parameter_name:
            p_boolean = []
            for flux_name in flux_parameter_list:
                boolean_existence = [True if j_parameter==parameter_k_flux else False
                                     for parameter_k_flux in flux_parameter_list[flux_name]]
                p_boolean.append(boolean_existence)
                if any(boolean_existence):
                    chosen_flux_name.append(flux_name)
                    chosen_flux_parameter_id.append([i for i, val in enumerate(boolean_existence) if val][0])
                    chosen_flux_orig_id.append(orig_id)
    else:
        for flux_name in flux_parameter_list:
            boolean_existence = [True if parameter_name[1] == parameter_k_flux else False
                                 for parameter_k_flux in flux_parameter_list[flux_name]]
            if any(boolean_existence):
                chosen_flux_name.append(flux_name)
                chosen_flux_parameter_id.append([i for i, val in enumerate(boolean_existence) if val])
                chosen_flux_orig_id.append(parameter_name[0])
    return chosen_flux_name, chosen_flux_parameter_id, chosen_flux_orig_id


def kotte_parameter_name(parameter_id):
    parameter_list = ['K1ac', 'K3fdp', 'L3fdp', 'K3pep',
                      'K2pep', 'vemax', 'Kefdp', 'ne', 'd',
                      'V4max', 'k1cat', 'V3max', 'V2max', 'ac']
    try:
        return [parameter_list[id] for id in parameter_id]
    except TypeError:
        return parameter_list[parameter_id]


def kotte_experiment_type_name(experiment_id):
    experiment_type_name_list = ['wildtype acetate', 'acetate perturbation', 'k1cat perturbation',
                                 'V3max perturbation', 'V2max perturbation']
    try:
        return [experiment_type_name_list[index] for index in experiment_id]
    except TypeError:
        return experiment_type_name_list[experiment_id]


def experiment_name(experiment_id, experiment_details):
    try:
        parameter_changed = kotte_parameter_name([int(experiment_details["indices"][i, 0]) for i in experiment_id])
        parameter_value = [experiment_details["indices"][i, 1] for i in experiment_id]
    except TypeError:
        parameter_changed = kotte_parameter_name(int(experiment_details["indices"][experiment_id, 0]))
        parameter_value = experiment_details["indices"][experiment_id, 1]
    experiment_name_list = ['{} changes {}'.format(j_p_change, j_p_value)
                            for j_p_change, j_p_value in zip(parameter_changed, parameter_value)]
    return experiment_name_list


def process_ident_data(ident_values, number_data):
    # create signed boolean array for identifiability
    signed_ident_values = np.sign(ident_values)
    ident_fun_val = []
    for id in range(0, number_data):
        ident_fun_val.append(signed_ident_values[id * 12:(id + 1) * 12, -1])
    p_list = [[p_id for p_id, val in enumerate(data_set) if val > 0] for data_set in ident_fun_val]
    p_list_boolean = [[True if parameter_id in list_1 else False for parameter_id in range(0, 12)] for list_1 in p_list]
    return p_list, np.array(p_list_boolean)


def boolean_ident_info(ident_values, number_of_parameters):
    signed_ident_values = np.sign(ident_values)
    p_list = [[p_id for p_id, val in enumerate(data_set) if val > 0] for data_set in signed_ident_values]
    p_list_boolean = [[True if parameter_id in list_1 else False for parameter_id in range(0, number_of_parameters)]
                      for list_1 in p_list]
    return p_list, np.array(p_list_boolean)


def create_data_for_file(ident_details, number_fluxes, fp_list, data_list):
    """create data for write/append all data to file"""
    number_data, p = ident_details["boolean"].shape
    write_2_file_data = []
    write_2_file_data.append(['Identifiable Perturbations'])
    flux_id, parameter_id, i = 1, 1, 0
    # perturbation_keys = parameter_list.keys()
    while i < p:
        write_2_file_data.append(['Flux {}, Parameter {}'.format(flux_id, parameter_id)])
        key_id = 'flux{}p{}'.format(flux_id, parameter_id)
        write_2_file_data.append(fp_list[key_id])
        if parameter_id < ident_details["parameters"][flux_id - 1]:
            parameter_id += 1
        else:
            if flux_id < number_fluxes:
                flux_id += 1
            else:
                flux_id = 1
            parameter_id = 1
        i += 1
    for j in range(0, number_data):
        write_2_file_data.append(['Parameters Identified by Data set {}'.format(j + 1)])
        key_id = 'dataset{}'.format(j + 1)
        write_2_file_data.append(data_list[key_id])
        # write_2_file_data.append(['Parameter Perturbations in Data Set {}'.format(j+1)])
        # write_2_file_data.append(perturbation_list[key_id])
    return write_2_file_data


def write_results_2_file(ident_details, number_fluxes, fp_list, data_list):
    write_2_file_data = create_data_for_file(ident_details, number_fluxes, fp_list, data_list)
    # write results to file
    path = "~" + "shyam" + r"\Documents\Courses\CHE1125Project\Results\ident\python2\kotte_ident_results.txt"
    fullpath = os.path.expanduser(path)
    with open(fullpath, 'a') as fh:
       writer = csv.writer(fh, delimiter="\t")
       [writer.writerow(r) for r in write_2_file_data]
       fh.close()
    return None


    # dataset dependent classification of parameters
    # data_list = data_based_processing(ident_details)

    # most useful dataset - based on number of parameter identified
    # max_data = get_most_useful_dataset(ident_details["boolean"])


    # choose additional data sets and consequently, experiments (if possible) to identify other parameters not
    # identified by chosen data set(s)
    # new_combos = calculate_experiment_combos(ident_details, experiment_details, perturbation_details, data_list)


def flux_ident_2_data_combination(all_data, flux_ids, choose=(), flux_choice=(), ident_fun_choice=()):
    """perform identifiability separately for each set of functions and generate separate identifiability info"""
    # 2 data combination ident list
    if flux_choice[0] == 1:
        flux_1 = flux_1_Vmax_ident
        flux_5 = flux_5_value1_ident
        flux_6 = flux_3_var1
    elif flux_choice[0] == 2:
        flux_1 = flux_1_kcat_ident
        flux_5 = flux_5_value2_ident
        flux_6 = flux_3_var2
    else:
        flux_1 = flux_1_ident_expression
        flux_5 = flux_5_ident_expression
        flux_6 = flux_3_var_1_and_2

    all_ident_fun_2_data = (flux_1, flux_2_ident_expression, flux_5, flux_6)

    # choose which flux functions to test identifiability for
    if ident_fun_choice:
        ident_fun_2_data = [all_ident_fun_2_data[j_ident_fun] for j_ident_fun in ident_fun_choice]
    else:
        ident_fun_2_data = all_ident_fun_2_data

    # run flux identifiability
    all_sample_all_fun_ident_info = multi_sample_ident_fun(ident_fun_2_data, all_data, choose, flux_ids, flux_choice)
    return all_sample_all_fun_ident_info


def flux_ident_3_data_combination(all_data, flux_ids, choose=(), flux_choice=(), ident_fun_choice=()):
    """perform identifiability separately for each set of functions and generate separate identifiability info"""
    # 3 data combination ident list
    if flux_choice[0] == 1:
        flux_3 = flux_3_value1_ident
    elif flux_choice[0] == 2:
        flux_3 = flux_3_value2_ident
    else:
        flux_3 = flux_3_ident_expression
    all_ident_fun_3_data = [flux_3]

    # choose flux functions to test identifiability for
    if ident_fun_choice:
        if len(ident_fun_choice) == len(all_ident_fun_3_data):
            ident_fun_3_data = all_ident_fun_3_data
        else:
            try:
                ident_fun_3_data = [all_ident_fun_3_data[j_ident_fun] for j_ident_fun in ident_fun_choice]
            except TypeError:
                ident_fun_3_data = all_ident_fun_3_data[ident_fun_choice]
    else:
        ident_fun_3_data = all_ident_fun_3_data
    all_sample_all_fun_ident_info = multi_sample_ident_fun(ident_fun_3_data, all_data, choose, flux_ids, flux_choice)
    return all_sample_all_fun_ident_info


def establish_kotte_flux_identifiability(all_data, choose=()):
    """call all identifiability evaluation funcs above and print numerical results"""
    # identifiability function list
    ident_fun_list = (flux_1_ident_expression, flux_2_ident_expression, flux_3_ident_expression)

    all_sample_ident_details = []
    for i_sample, sample_data in enumerate(all_data):
        if choose:
            try:
                chosen_values = list(sample_data["values"][:choose, :])
            except TypeError:
                iter_chosen_value = []
                for indexes in range(0, len(choose)):
                    iter_chosen_value.append(list(sample_data["values"][indexes, :]))
                chosen_values = iter_chosen_value[:]
        else:
            chosen_values = list(sample_data["values"][:, :])
        # chosen_details = all_data["details"][0:choose, :]

        # number_fluxes = len(ident_fun_list)
        number_data = len(chosen_values)
        ident_values, parameters_per_flux, all_fun_array_list = get_ident_value(ident_fun_list, chosen_values, choose)

        # process identifiability data
        _, plist_boolean = process_ident_data(ident_values, number_data)
        ident_details = {"boolean": plist_boolean, "values": ident_values, "parameters": parameters_per_flux}

        # flux_name, parameter_id = return_flux_based_id([5])
        # parameter_name = ident_parameter_name(parameter_id, flux_name)
        all_sample_ident_details.append(ident_details)

    return all_sample_ident_details


def select_experiment_for_dataset(experiment_id):
    return None

def loop_through_experiments(experiments_per_dataset, total_experiments):
    chosen_ids = [0]
    for experiment in experiments_per_dataset:

        pass
    experimental_data = []
    return experimental_data
