import numpy as np
import sympy as sym
import itertools as it
import os.path
import csv
# import scipy.linalg

# K1ac, K3fdp, L3fdp, K3pep, K2pep, vemax, Kefdp, ne, d, V4max, k1cat, V3max, V2max, ac
def_par_val = np.array([.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, .1])


def kotte_ck_flux(y, p=def_par_val):
    """calculate flux using convenience kinetics"""

    K1ac, K3fdp, L3fdp, K3pep, K2pep, vemax, Kefdp, ne, d, V4max, k1cat, V3max, V2max, ac = p

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


def kotte_flux(y, p=def_par_val):
    """function doc_string"""

    K1ac, K3fdp, L3fdp, K3pep, K2pep, vemax, Kefdp, ne, d, V4max, k1cat, V3max, V2max, ac = p

    flux_1 = k1cat*y[2]*ac/(ac+K1ac)
    flux_2 = vemax*(1-1/(1+(Kefdp/y[1])**ne))
    fdp_sat = 1 + y[1]/K3fdp
    pep_sat = 1 + y[0]/K3pep
    flux_3 = V3max*(fdp_sat-1)*(fdp_sat)**3/(fdp_sat**4+L3fdp*(pep_sat**(-4)))
    flux_4 = V2max*y[0]/(y[0]+K2pep)
    flux_5 = V4max*y[0]
    flux_6 = d*y[2]
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


def define_sym_variables():
    """define all required symbolic variables for sympy expressions"""
    ac1, ac2, ac3, x11, x12, x13, x21, x22, x23, x31, x32, x33, \
    v31, v32, v33, v11, v12, v13, v21, v22, v23, v41, v42, v43 = \
        sym.symbols('ac1, ac2, ac3, x11, x12, x13, x21, x22, x23, x31, x32, x33,'
                ' v31, v32, v33, v11, v12, v13, v21, v22, v23, v41, v42, v43', positive=True)
    variables = [ac1, x11, x21, x31, v11, v21, v31, v41,
                 ac2, x12, x22, x32, v12, v22, v32, v42,
                 ac3, x13, x23, x33, v13, v23, v33, v43]
    return variables, ac1, ac2, ac3, x11, x12, x13, x21, x22, x23, x31, x32, x33, \
           v31, v32, v33, v11, v12, v13, v21, v22, v23, v41, v42, v43


def flux_1_ident_expression(experimental_data):
    """symbolic and lambdify expression for flux 1 denominator from mathematica"""
    # define symbols and variables (obtained through experimental data
    variables, ac1, ac2, _,\
        _, _, _, _, _, _, x31, x32, _,\
        _, _, _, v11, v12, _, \
        _, _, _, _, _, _ = define_sym_variables()

    # get variable values (w/o sympy directly from experimental data)
    ac1, x11, x21, x31, v11, v21, v31, v41, \
    ac2, x12, x22, x32, v12, v22, v32, v42, \
    ac3, x13, x23, x33, v13, v23, v33, v43 = list(experimental_data)
    #ac = experimental_data[range(0, 24, 8)]
    #x_ranges = [range(i, 24, 8) for i in range(1, 4)]
    #all_x_ranges = [list(np.array(x_ranges)[:, i]) for i in range(0, 3)]
    #x = [experimental_data[all_x_ranges[i]] for i in range(0, 3)]
    #v_ranges = [range(i, 24, 8) for i in range(4, 8)]
    #all_v_ranges = [list(np.array(v_ranges)[:, j]) for j in range(0, 3)]
    #v = [experimental_data[all_v_ranges[j]] for j in range(0, 3)]

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
    # define symbols and variables (obtained through experimental data
    variables, _, _, _, \
    _, _, _, x21, x22, _, _, _, _, \
    _, _, _, _, _, _, \
    v21, v22, _, _, _, _ = define_sym_variables()

    # symbolic expression for v2
    # V2max
    v2max_numerator = -v21 * v22 * x21 + v21 * v22 * x22
    v2max_denominator = -(v22 * x21 - v21 * x22)
    # K2pep
    k2pep_numerator = x21 * (v21 * x22 - v22 * x22)
    k2pep_denominator = v22 * x21 - v21 * x22
    v2max_nr_expr = sym.lambdify([variables], v2max_numerator, "numpy")
    v2max_dr_expr = sym.lambdify([variables], v2max_denominator, "numpy")
    v2max_expr = sym.lambdify([variables], v2max_numerator/v2max_denominator, "numpy")
    k2pep_nr_expr = sym.lambdify([variables], k2pep_numerator, "numpy")
    k2pep_dr_expr = sym.lambdify([variables], k2pep_denominator, "numpy")
    k2pep_expr = sym.lambdify([variables], k2pep_numerator/k2pep_denominator, "numpy")
    v2max_numerator_value = v2max_nr_expr(experimental_data)
    v2max_denominator_value = v2max_dr_expr(experimental_data)
    v2max_value = v2max_expr(experimental_data)
    k2pep_numerator_value = k2pep_nr_expr(experimental_data)
    k2pep_denominator_value = k2pep_dr_expr(experimental_data)
    k2pep_value = k2pep_expr(experimental_data)

    return [v2max_numerator_value, v2max_denominator_value, v2max_value], \
           [k2pep_numerator_value, k2pep_denominator_value, k2pep_value]


def flux_3_ident_expression(experimental_data):
    """symbolic and lambdify expression for flux 3 denominator from mathematica"""
    # define symbols and variables (obtained through experimental data
    variables, _, _, _, \
    x11, x12, x13, x21, x22, x23, _, _, _, \
    v31, v32, v33, _, _, _, \
    _, _, _, _, _, _ = define_sym_variables()

    # symbolic expression for v3
    # V3max
    v3max_numerator_1 = -v31*v32*v33*x11*x12*x21 + v31*v32*v33*x11*x13*x21 + v31*v32*v33*x11*x12*x22 - \
                      v31*v32*v33*x12*x13*x22 - v31*v32*v33*x11*x13*x23 + v31*v32*v33*x12*x13*x23 - \
                      (v31*v32*v33*x12*x21*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                            v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                            v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                            v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                            v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                            v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 -
                                            sym.sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
                                                  v32*v33*x11*x13*x21*x22 + v31*v33*x12*x13*x21*x22 +
                                                  v32*v33*x11*x12*x21*x23 - v31*v32*x11*x13*x21*x23 +
                                                  v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x21*x23 -
                                                  v31*v33*x11*x12*x22*x23 + v31*v32*x11*x13*x22*x23 +
                                                  v31*v32*x12*x13*x22*x23 - v31*v33*x12*x13*x22*x23)**2 -
                                                 4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 -
                                                    v31*v32*x13*x21*x23 - v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                                                 (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 -
                                                  v31*v32*x11*x12*x13*x21*x23 + v32*v33*x11*x12*x13*x21*x23 +
                                                  v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))))/\
                      (2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                          v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)) + \
                      (v31*v32*v33*x13*x21*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                            v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                            v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                            v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                            v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                            v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 -
                                            sym.sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
                                                  v32*v33*x11*x13*x21*x22 + v31*v33*x12*x13*x21*x22 +
                                                  v32*v33*x11*x12*x21*x23 - v31*v32*x11*x13*x21*x23 +
                                                  v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x21*x23 -
                                                  v31*v33*x11*x12*x22*x23 + v31*v32*x11*x13*x22*x23 +
                                                  v31*v32*x12*x13*x22*x23 - v31*v33*x12*x13*x22*x23)**2 -
                                                 4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 -
                                                    v31*v32*x13*x21*x23 - v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                                                 (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 -
                                                  v31*v32*x11*x12*x13*x21*x23 + v32*v33*x11*x12*x13*x21*x23 +
                                                  v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))))/\
                      (2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                          v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)) + \
                      (v31*v32*v33*x11*x22*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                            v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                            v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                            v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                            v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                            v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 -
                                            sym.sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
                                                  v32*v33*x11*x13*x21*x22 + v31*v33*x12*x13*x21*x22 +
                                                  v32*v33*x11*x12*x21*x23 - v31*v32*x11*x13*x21*x23 +
                                                  v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x21*x23 -
                                                  v31*v33*x11*x12*x22*x23 + v31*v32*x11*x13*x22*x23 +
                                                  v31*v32*x12*x13*x22*x23 - v31*v33*x12*x13*x22*x23)**2 -
                                                 4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 -
                                                    v31*v32*x13*x21*x23 - v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                                                 (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 -
                                                  v31*v32*x11*x12*x13*x21*x23 + v32*v33*x11*x12*x13*x21*x23 +
                                                  v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))))/\
                      (2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                          v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)) - \
                      (v31*v32*v33*x13*x22*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                            v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                            v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                            v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                            v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                            v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 -
                                            sym.sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
                                                  v32*v33*x11*x13*x21*x22 + v31*v33*x12*x13*x21*x22 +
                                                  v32*v33*x11*x12*x21*x23 - v31*v32*x11*x13*x21*x23 +
                                                  v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x21*x23 -
                                                  v31*v33*x11*x12*x22*x23 + v31*v32*x11*x13*x22*x23 +
                                                  v31*v32*x12*x13*x22*x23 - v31*v33*x12*x13*x22*x23)**2 -
                                                 4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 -
                                                    v31*v32*x13*x21*x23 - v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                                                 (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 -
                                                  v31*v32*x11*x12*x13*x21*x23 + v32*v33*x11*x12*x13*x21*x23 +
                                                  v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))))/\
                      (2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                          v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)) - \
                      (v31*v32*v33*x11*x23*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                            v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                            v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                            v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                            v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                            v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 -
                                            sym.sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
                                                  v32*v33*x11*x13*x21*x22 + v31*v33*x12*x13*x21*x22 +
                                                  v32*v33*x11*x12*x21*x23 - v31*v32*x11*x13*x21*x23 +
                                                  v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x21*x23 -
                                                  v31*v33*x11*x12*x22*x23 + v31*v32*x11*x13*x22*x23 +
                                                  v31*v32*x12*x13*x22*x23 - v31*v33*x12*x13*x22*x23)**2 -
                                                 4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 -
                                                    v31*v32*x13*x21*x23 - v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                                                 (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 -
                                                  v31*v32*x11*x12*x13*x21*x23 + v32*v33*x11*x12*x13*x21*x23 +
                                                  v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))))/\
                      (2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                          v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)) + \
                      (v31*v32*v33*x12*x23*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                            v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                            v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                            v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                            v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                            v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 -
                                            sym.sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
                                                  v32*v33*x11*x13*x21*x22 + v31*v33*x12*x13*x21*x22 +
                                                  v32*v33*x11*x12*x21*x23 - v31*v32*x11*x13*x21*x23 +
                                                  v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x21*x23 -
                                                  v31*v33*x11*x12*x22*x23 + v31*v32*x11*x13*x22*x23 +
                                                  v31*v32*x12*x13*x22*x23 - v31*v33*x12*x13*x22*x23)**2 -
                                                 4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 -
                                                    v31*v32*x13*x21*x23 - v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                                                 (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 -
                                                  v31*v32*x11*x12*x13*x21*x23 + v32*v33*x11*x12*x13*x21*x23 +
                                                  v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))))/\
                      (2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                          v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23))
    v3max_denominator_1 = -v32*v33*x11*x12*x21 + v32*v33*x11*x13*x21 + v31*v33*x11*x12*x22 -  \
                        v31*v33*x12*x13*x22 - v31*v32*x11*x13*x23 + v31*v32*x12*x13*x23
    v3max_nr_1_expr = sym.lambdify([variables], v3max_numerator_1, "numpy")
    v3max_dr_1_expr = sym.lambdify([variables], v3max_denominator_1, "numpy")
    v3max_1_expr = sym.lambdify([variables], v3max_numerator_1 / v3max_denominator_1, "numpy")
    v3max_nr_1_value = v3max_nr_1_expr(experimental_data)
    v3max_dr_1_value = v3max_dr_1_expr(experimental_data)
    v3max_1_value = v3max_1_expr(experimental_data)

    # v3max = second solution
    v3max_numerator_2 = -v31*v32*v33*x11*x12*x21 + v31*v32*v33*x11*x13*x21 + v31*v32*v33*x11*x12*x22 - \
                        v31*v32*v33*x12*x13*x22 - v31*v32*v33*x11*x13*x23 + v31*v32*v33*x12*x13*x23 - \
                        (v31*v32*v33*x12*x21*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                              v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                              v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                              v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                              v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                              v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 +
                                              sym.sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
                                                    v32*v33*x11*x13*x21*x22 + v31*v33*x12*x13*x21*x22 +
                                                    v32*v33*x11*x12*x21*x23 - v31*v32*x11*x13*x21*x23 +
                                                    v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x21*x23 -
                                                    v31*v33*x11*x12*x22*x23 + v31*v32*x11*x13*x22*x23 +
                                                    v31*v32*x12*x13*x22*x23 - v31*v33*x12*x13*x22*x23)**2 -
                                                   4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 -
                                                      v31*v32*x13*x21*x23 - v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                                                   (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 -
                                                    v31*v32*x11*x12*x13*x21*x23 + v32*v33*x11*x12*x13*x21*x23 +
                                                    v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))))/\
                        (2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                            v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)) + \
                        (v31*v32*v33*x13*x21*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                              v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                              v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                              v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                              v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                              v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 +
                                              sym.sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
                                                    v32*v33*x11*x13*x21*x22 + v31*v33*x12*x13*x21*x22 +
                                                    v32*v33*x11*x12*x21*x23 - v31*v32*x11*x13*x21*x23 +
                                                    v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x21*x23 -
                                                    v31*v33*x11*x12*x22*x23 + v31*v32*x11*x13*x22*x23 +
                                                    v31*v32*x12*x13*x22*x23 - v31*v33*x12*x13*x22*x23)**2 -
                                                   4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 -
                                                      v31*v32*x13*x21*x23 - v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                                                   (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 -
                                                    v31*v32*x11*x12*x13*x21*x23 + v32*v33*x11*x12*x13*x21*x23 +
                                                    v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))))/\
                        (2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                            v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)) + \
                        (v31*v32*v33*x11*x22*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                              v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                              v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                              v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                              v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                              v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 +
                                              sym.sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
                                                    v32*v33*x11*x13*x21*x22 + v31*v33*x12*x13*x21*x22 +
                                                    v32*v33*x11*x12*x21*x23 - v31*v32*x11*x13*x21*x23 +
                                                    v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x21*x23 -
                                                    v31*v33*x11*x12*x22*x23 + v31*v32*x11*x13*x22*x23 +
                                                    v31*v32*x12*x13*x22*x23 - v31*v33*x12*x13*x22*x23)**2 -
                                                   4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 -
                                                      v31*v32*x13*x21*x23 - v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                                                   (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 -
                                                    v31*v32*x11*x12*x13*x21*x23 + v32*v33*x11*x12*x13*x21*x23 +
                                                    v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))))/\
                        (2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                            v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)) - \
                        (v31*v32*v33*x13*x22*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                              v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                              v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                              v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                              v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                              v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 +
                                              sym.sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
                                                    v32*v33*x11*x13*x21*x22 + v31*v33*x12*x13*x21*x22 +
                                                    v32*v33*x11*x12*x21*x23 - v31*v32*x11*x13*x21*x23 +
                                                    v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x21*x23 -
                                                    v31*v33*x11*x12*x22*x23 + v31*v32*x11*x13*x22*x23 +
                                                    v31*v32*x12*x13*x22*x23 - v31*v33*x12*x13*x22*x23)**2 -
                                                   4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 -
                                                      v31*v32*x13*x21*x23 - v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                                                   (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 -
                                                    v31*v32*x11*x12*x13*x21*x23 + v32*v33*x11*x12*x13*x21*x23 +
                                                    v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))))/\
                        (2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                            v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)) - \
                        (v31*v32*v33*x11*x23*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                              v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                              v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                              v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                              v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                              v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 +
                                              sym.sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
                                                    v32*v33*x11*x13*x21*x22 + v31*v33*x12*x13*x21*x22 +
                                                    v32*v33*x11*x12*x21*x23 - v31*v32*x11*x13*x21*x23 +
                                                    v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x21*x23 -
                                                    v31*v33*x11*x12*x22*x23 + v31*v32*x11*x13*x22*x23 +
                                                    v31*v32*x12*x13*x22*x23 - v31*v33*x12*x13*x22*x23)**2 -
                                                   4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 -
                                                      v31*v32*x13*x21*x23 - v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                                                   (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 -
                                                    v31*v32*x11*x12*x13*x21*x23 + v32*v33*x11*x12*x13*x21*x23 +
                                                    v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))))/\
                        (2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                            v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)) + \
                        (v31*v32*v33*x12*x23*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                              v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                              v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                              v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                              v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                              v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 +
                                              sym.sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
                                                    v32*v33*x11*x13*x21*x22 + v31*v33*x12*x13*x21*x22 +
                                                    v32*v33*x11*x12*x21*x23 - v31*v32*x11*x13*x21*x23 +
                                                    v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x21*x23 -
                                                    v31*v33*x11*x12*x22*x23 + v31*v32*x11*x13*x22*x23 +
                                                    v31*v32*x12*x13*x22*x23 - v31*v33*x12*x13*x22*x23)**2 -
                                                   4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 -
                                                      v31*v32*x13*x21*x23 - v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                                                   (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 -
                                                    v31*v32*x11*x12*x13*x21*x23 + v32*v33*x11*x12*x13*x21*x23 +
                                                    v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))))/\
                        (2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                            v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23))
    v3max_denominator_2 = -v32*v33*x11*x12*x21 + v32*v33*x11*x13*x21 + v31*v33*x11*x12*x22 -  \
                          v31*v33*x12*x13*x22 - v31*v32*x11*x13*x23 + v31*v32*x12*x13*x23
    v3max_nr_2_expr = sym.lambdify([variables], v3max_numerator_2, "numpy")
    v3max_dr_2_expr = sym.lambdify([variables], v3max_denominator_2, "numpy")
    v3max_2_expr = sym.lambdify([variables], v3max_numerator_2/v3max_denominator_2, "numpy")
    v3max_nr_2_value = v3max_nr_2_expr(experimental_data)
    v3max_dr_2_value = v3max_dr_2_expr(experimental_data)
    v3max_2_value = v3max_2_expr(experimental_data)

    # K3fdp
    k3fdp_numerator_1 = -v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 + v31*v32*x11*x13*x21*x23 - \
                        v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 + \
                        (v32*v33*x11*x21*x22*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                              v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                              v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                              v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                              v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                              v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 -
                                              sym.sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
                                                    v32*v33*x11*x13*x21*x22 + v31*v33*x12*x13*x21*x22 +
                                                    v32*v33*x11*x12*x21*x23 - v31*v32*x11*x13*x21*x23 +
                                                    v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x21*x23 -
                                                    v31*v33*x11*x12*x22*x23 + v31*v32*x11*x13*x22*x23 +
                                                    v31*v32*x12*x13*x22*x23 - v31*v33*x12*x13*x22*x23)**2 -
                                                   4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 +
                                                      v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                                                      v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                                                   (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 -
                                                    v31*v32*x11*x12*x13*x21*x23 + v32*v33*x11*x12*x13*x21*x23 +
                                                    v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))))/\
                        (2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                            v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)) - \
                        (v31*v33*x12*x21*x22*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                              v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                              v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                              v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                              v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                              v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 -
                                              sym.sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
                                                    v32*v33*x11*x13*x21*x22 + v31*v33*x12*x13*x21*x22 +
                                                    v32*v33*x11*x12*x21*x23 - v31*v32*x11*x13*x21*x23 +
                                                    v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x21*x23 -
                                                    v31*v33*x11*x12*x22*x23 + v31*v32*x11*x13*x22*x23 +
                                                    v31*v32*x12*x13*x22*x23 - v31*v33*x12*x13*x22*x23)**2 -
                                                   4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 +
                                                      v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                                                      v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                                                   (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 -
                                                    v31*v32*x11*x12*x13*x21*x23 + v32*v33*x11*x12*x13*x21*x23 +
                                                    v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))))/ \
                        (2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                            v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)) - \
                        (v32*v33*x11*x21*x23*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                              v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                              v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                              v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                              v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                              v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 -
                                              sym.sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
                                                    v32*v33*x11*x13*x21*x22 + v31*v33*x12*x13*x21*x22 +
                                                    v32*v33*x11*x12*x21*x23 - v31*v32*x11*x13*x21*x23 +
                                                    v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x21*x23 -
                                                    v31*v33*x11*x12*x22*x23 + v31*v32*x11*x13*x22*x23 +
                                                    v31*v32*x12*x13*x22*x23 - v31*v33*x12*x13*x22*x23)**2 -
                                                   4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 +
                                                      v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                                                      v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                                                   (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 -
                                                    v31*v32*x11*x12*x13*x21*x23 + v32*v33*x11*x12*x13*x21*x23 +
                                                    v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))))/ \
                        (2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 -
                            v31*v32*x13*x21*x23 - v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)) + \
                        (v31*v32*x13*x21*x23*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                              v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                              v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                              v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                              v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                              v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 -
                                              sym.sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
                                                    v32*v33*x11*x13*x21*x22 + v31*v33*x12*x13*x21*x22 +
                                                    v32*v33*x11*x12*x21*x23 - v31*v32*x11*x13*x21*x23 +
                                                    v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x21*x23 -
                                                    v31*v33*x11*x12*x22*x23 + v31*v32*x11*x13*x22*x23 +
                                                    v31*v32*x12*x13*x22*x23 - v31*v33*x12*x13*x22*x23)**2 -
                                                   4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 +
                                                      v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                                                      v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                                                   (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 -
                                                    v31*v32*x11*x12*x13*x21*x23 + v32*v33*x11*x12*x13*x21*x23 +
                                                    v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))))/ \
                        (2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                            v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)) + \
                        (v31*v33*x12*x22*x23*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                              v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                              v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                              v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                              v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                              v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 -
                                              sym.sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
                                                    v32*v33*x11*x13*x21*x22 + v31*v33*x12*x13*x21*x22 +
                                                    v32*v33*x11*x12*x21*x23 - v31*v32*x11*x13*x21*x23 +
                                                    v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x21*x23 -
                                                    v31*v33*x11*x12*x22*x23 + v31*v32*x11*x13*x22*x23 +
                                                    v31*v32*x12*x13*x22*x23 - v31*v33*x12*x13*x22*x23)**2 -
                                                   4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 +
                                                      v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                                                      v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                                                   (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 -
                                                    v31*v32*x11*x12*x13*x21*x23 + v32*v33*x11*x12*x13*x21*x23 +
                                                    v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))))/ \
                        (2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                            v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)) - \
                        (v31*v32*x13*x22*x23*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                              v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                              v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                              v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                              v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                              v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 -
                                              sym.sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
                                                    v32*v33*x11*x13*x21*x22 + v31*v33*x12*x13*x21*x22 +
                                                    v32*v33*x11*x12*x21*x23 - v31*v32*x11*x13*x21*x23 +
                                                    v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x21*x23 -
                                                    v31*v33*x11*x12*x22*x23 + v31*v32*x11*x13*x22*x23 +
                                                    v31*v32*x12*x13*x22*x23 - v31*v33*x12*x13*x22*x23)**2 -
                                                   4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 +
                                                      v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                                                      v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                                                   (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 -
                                                    v31*v32*x11*x12*x13*x21*x23 + v32*v33*x11*x12*x13*x21*x23 +
                                                    v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))))/ \
                        (2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                            v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23))
    k3fdp_denominator_1 = -v32*v33*x11*x12*x21 + v32*v33*x11*x13*x21 + v31*v33*x11*x12*x22 -  \
                          v31*v33*x12*x13*x22 - v31*v32*x11*x13*x23 + v31*v32*x12*x13*x23
    k3fdp_nr_1_expr = sym.lambdify([variables], k3fdp_numerator_1, "numpy")
    k3fdp_dr_1_expr = sym.lambdify([variables], k3fdp_denominator_1, "numpy")
    k3fdp_1_expr = sym.lambdify([variables], k3fdp_numerator_1 / k3fdp_denominator_1, "numpy")
    k3fdp_nr_1_value = k3fdp_nr_1_expr(experimental_data)
    k3fdp_dr_1_value = k3fdp_dr_1_expr(experimental_data)
    k3fdp_1_value = k3fdp_1_expr(experimental_data)

    # K3fdp 2
    k3fdp_numerator_2 = -v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 + v31*v32*x11*x13*x21*x23 - \
                        v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 + \
                        (v32*v33*x11*x21*x22*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                              v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                              v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                              v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                              v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                              v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 +
                                              sym.sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
                                                    v32*v33*x11*x13*x21*x22 + v31*v33*x12*x13*x21*x22 +
                                                    v32*v33*x11*x12*x21*x23 - v31*v32*x11*x13*x21*x23 +
                                                    v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x21*x23 -
                                                    v31*v33*x11*x12*x22*x23 + v31*v32*x11*x13*x22*x23 +
                                                    v31*v32*x12*x13*x22*x23 - v31*v33*x12*x13*x22*x23)**2 -
                                                   4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 -
                                                      v31*v32*x13*x21*x23 - v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                                                   (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 -
                                                    v31*v32*x11*x12*x13*x21*x23 + v32*v33*x11*x12*x13*x21*x23 +
                                                    v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))))/\
                        (2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                            v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)) - \
                        (v31*v33*x12*x21*x22*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                              v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                              v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                              v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                              v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                              v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 +
                                              sym.sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
                                                    v32*v33*x11*x13*x21*x22 + v31*v33*x12*x13*x21*x22 +
                                                    v32*v33*x11*x12*x21*x23 - v31*v32*x11*x13*x21*x23 +
                                                    v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x21*x23 -
                                                    v31*v33*x11*x12*x22*x23 + v31*v32*x11*x13*x22*x23 +
                                                    v31*v32*x12*x13*x22*x23 - v31*v33*x12*x13*x22*x23)**2 -
                                                   4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 -
                                                      v31*v32*x13*x21*x23 - v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                                                   (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 -
                                                    v31*v32*x11*x12*x13*x21*x23 + v32*v33*x11*x12*x13*x21*x23 +
                                                    v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))))/\
                        (2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                            v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)) - \
                        (v32*v33*x11*x21*x23*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                              v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                              v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                              v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                              v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                              v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 +
                                              sym.sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
                                                    v32*v33*x11*x13*x21*x22 + v31*v33*x12*x13*x21*x22 +
                                                    v32*v33*x11*x12*x21*x23 - v31*v32*x11*x13*x21*x23 +
                                                    v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x21*x23 -
                                                    v31*v33*x11*x12*x22*x23 + v31*v32*x11*x13*x22*x23 +
                                                    v31*v32*x12*x13*x22*x23 - v31*v33*x12*x13*x22*x23)**2 -
                                                   4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 -
                                                      v31*v32*x13*x21*x23 - v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                                                   (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 -
                                                    v31*v32*x11*x12*x13*x21*x23 + v32*v33*x11*x12*x13*x21*x23 +
                                                    v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))))/\
                        (2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                            v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)) + \
                        (v31*v32*x13*x21*x23*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                              v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                              v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                              v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                              v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                              v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 +
                                              sym.sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
                                                    v32*v33*x11*x13*x21*x22 + v31*v33*x12*x13*x21*x22 +
                                                    v32*v33*x11*x12*x21*x23 - v31*v32*x11*x13*x21*x23 +
                                                    v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x21*x23 -
                                                    v31*v33*x11*x12*x22*x23 + v31*v32*x11*x13*x22*x23 +
                                                    v31*v32*x12*x13*x22*x23 - v31*v33*x12*x13*x22*x23)**2 -
                                                   4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 -
                                                      v31*v32*x13*x21*x23 - v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                                                   (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 -
                                                    v31*v32*x11*x12*x13*x21*x23 + v32*v33*x11*x12*x13*x21*x23 +
                                                    v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))))/\
                        (2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                            v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)) + \
                        (v31*v33*x12*x22*x23*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                              v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                              v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                              v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                              v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                              v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 +
                                              sym.sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
                                                    v32*v33*x11*x13*x21*x22 + v31*v33*x12*x13*x21*x22 +
                                                    v32*v33*x11*x12*x21*x23 - v31*v32*x11*x13*x21*x23 +
                                                    v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x21*x23 -
                                                    v31*v33*x11*x12*x22*x23 + v31*v32*x11*x13*x22*x23 +
                                                    v31*v32*x12*x13*x22*x23 - v31*v33*x12*x13*x22*x23)**2 -
                                                   4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 -
                                                      v31*v32*x13*x21*x23 - v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                                                   (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 -
                                                    v31*v32*x11*x12*x13*x21*x23 + v32*v33*x11*x12*x13*x21*x23 +
                                                    v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))))/\
                        (2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                            v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)) - \
                        (v31*v32*x13*x22*x23*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                              v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                              v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                              v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                              v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                              v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 +
                                              sym.sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
                                                    v32*v33*x11*x13*x21*x22 + v31*v33*x12*x13*x21*x22 +
                                                    v32*v33*x11*x12*x21*x23 - v31*v32*x11*x13*x21*x23 +
                                                    v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x21*x23 -
                                                    v31*v33*x11*x12*x22*x23 + v31*v32*x11*x13*x22*x23 +
                                                    v31*v32*x12*x13*x22*x23 - v31*v33*x12*x13*x22*x23)**2 -
                                                   4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 -
                                                      v31*v32*x13*x21*x23 - v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                                                   (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 -
                                                    v31*v32*x11*x12*x13*x21*x23 + v32*v33*x11*x12*x13*x21*x23 +
                                                    v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))))/\
                        (2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                            v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23))
    k3fdp_denominator_2 = -v32*v33*x11*x12*x21 + v32*v33*x11*x13*x21 + v31*v33*x11*x12*x22 -  \
                          v31*v33*x12*x13*x22 - v31*v32*x11*x13*x23 + v31*v32*x12*x13*x23
    k3fdp_nr_2_expr = sym.lambdify([variables], k3fdp_numerator_2, "numpy")
    k3fdp_dr_2_expr = sym.lambdify([variables], k3fdp_denominator_2, "numpy")
    k3fdp_2_expr = sym.lambdify([variables], k3fdp_numerator_2 / k3fdp_denominator_2, "numpy")
    k3fdp_nr_2_value = k3fdp_nr_2_expr(experimental_data)
    k3fdp_dr_2_value = k3fdp_dr_2_expr(experimental_data)
    k3fdp_2_value = k3fdp_2_expr(experimental_data)

    # K3pep
    k3pep_numerator_1 = -v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 + v32*v33*x11*x13*x21*x22 - \
                        v31*v33*x12*x13*x21*x22 - v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 - \
                        v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 + v31*v33*x11*x12*x22*x23 - \
                        v31*v32*x11*x13*x22*x23 - v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 - \
                        sym.sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 - v32*v33*x11*x13*x21*x22 +
                              v31*v33*x12*x13*x21*x22 + v32*v33*x11*x12*x21*x23 - v31*v32*x11*x13*x21*x23 +
                              v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x21*x23 - v31*v33*x11*x12*x22*x23 +
                              v31*v32*x11*x13*x22*x23 + v31*v32*x12*x13*x22*x23 - v31*v33*x12*x13*x22*x23)**2 -
                             4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                                v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                             (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 - v31*v32*x11*x12*x13*x21*x23 +
                              v32*v33*x11*x12*x13*x21*x23 + v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))
    k3pep_denominator_1 = 2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 -
                             v31*v32*x13*x21*x23 - v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)
    k3pep_nr_1_expr = sym.lambdify([variables], k3pep_numerator_1, "numpy")
    k3pep_dr_1_expr = sym.lambdify([variables], k3pep_denominator_1, "numpy")
    k3pep_1_expr = sym.lambdify([variables], k3pep_numerator_1 / k3pep_denominator_1, "numpy")
    k3pep_nr_1_value = k3pep_nr_1_expr(experimental_data)
    k3pep_dr_1_value = k3pep_dr_1_expr(experimental_data)
    k3pep_1_value = k3pep_1_expr(experimental_data)

    k3pep_numerator_2 = -v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 + v32*v33*x11*x13*x21*x22 - \
                        v31*v33*x12*x13*x21*x22 - v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 - \
                        v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 + v31*v33*x11*x12*x22*x23 - \
                        v31*v32*x11*x13*x22*x23 - v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 + \
                        sym.sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 - v32*v33*x11*x13*x21*x22 +
                              v31*v33*x12*x13*x21*x22 + v32*v33*x11*x12*x21*x23 - v31*v32*x11*x13*x21*x23 +
                              v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x21*x23 - v31*v33*x11*x12*x22*x23 +
                              v31*v32*x11*x13*x22*x23 + v31*v32*x12*x13*x22*x23 - v31*v33*x12*x13*x22*x23)**2 -
                             4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                                v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                             (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 - v31*v32*x11*x12*x13*x21*x23 +
                              v32*v33*x11*x12*x13*x21*x23 + v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))
    k3pep_denominator_2 = 2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 -
                             v31*v32*x13*x21*x23 - v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)
    k3pep_nr_2_expr = sym.lambdify([variables], k3pep_numerator_2, "numpy")
    k3pep_dr_2_expr = sym.lambdify([variables], k3pep_denominator_2, "numpy")
    k3pep_2_expr = sym.lambdify([variables], k3pep_numerator_2 / k3pep_denominator_2, "numpy")
    k3pep_nr_2_value = k3pep_nr_2_expr(experimental_data)
    k3pep_dr_2_value = k3pep_dr_2_expr(experimental_data)
    k3pep_2_value = k3pep_2_expr(experimental_data)

    return [v3max_nr_1_value, v3max_dr_1_value, v3max_1_value], \
           [k3fdp_nr_1_value, k3fdp_dr_1_value, k3fdp_1_value], \
           [k3pep_nr_1_value, k3pep_dr_1_value, k3pep_1_value], \
           [v3max_nr_2_value, v3max_dr_2_value, v3max_2_value], \
           [k3fdp_nr_2_value, k3fdp_dr_2_value, k3fdp_2_value], \
           [k3pep_nr_2_value, k3pep_dr_2_value, k3pep_2_value]


def truncate_values(f, n=3):
    """truncates floats to n specified values after the decimal"""
    if not np.isnan(f):
        s = '{}'.format(f)  # convert float to string
        i, p, d = s.partition('.')
        return float('.'.join([i, (d+'0'*n)[:n]]))
    else:
        return f


def call_truncate_method(ident_value_list, parameter_count, expression_count=3):
    flux_ident_value = np.zeros((parameter_count, expression_count))
    for i, j in enumerate(ident_value_list):
        trunc_value = map(truncate_values, j)
        # trunc_value = map(float, trunc_value)
        flux_ident_value[i, :] = np.array(trunc_value)
    return flux_ident_value


def run_flux_ident(ident_function_list, data):
    number_of_parameters = 0
    number_of_parameters_per_flux = [None] * len(ident_function_list)
    ident_value_list = []
    iterator = 0
    for func in ident_function_list:
        ident_value = func(data)
        ident_value_list.append(ident_value)
        number_of_expressions = len(ident_value[0])
        number_of_parameters_per_flux[iterator] = len(ident_value)
        iterator += 1
        number_of_parameters += len(ident_value)

    ident_value_array = np.zeros((number_of_parameters, number_of_expressions))
    irow = 0
    for iflux in ident_value_list:
        truncated_ident_value = call_truncate_method(iflux, len(iflux))
        nrows, ncolumns = np.shape(truncated_ident_value)
        ident_value_array[irow:(irow + nrows), :] = truncated_ident_value
        irow += nrows
    return ident_value_array, number_of_parameters_per_flux, ncolumns


def get_ident_value(ident_function_list, experimental_data_list):
    all_data_ident_lists = []
    number_of_parameters_per_flux = [0]*len(ident_function_list)
    number_of_expressions_per_parameter = 0
    for index, data_set in enumerate(experimental_data_list):
        print('Identifiability for Dataset {} of {}\n'.format(index + 1, len(experimental_data_list)))
        identifiability_values, number_of_parameters_per_flux, number_of_expressions_per_parameter = \
            run_flux_ident(ident_function_list, data_set)
        all_data_ident_lists.append(identifiability_values)

    number_of_data_sets = len(experimental_data_list)
    total_parameters_identified = sum(number_of_parameters_per_flux)
    all_identifiability_values = \
        np.zeros((number_of_data_sets * total_parameters_identified, number_of_expressions_per_parameter))
    array_index = 0
    for idata in all_data_ident_lists:
        all_identifiability_values[array_index:(array_index + total_parameters_identified), :] = idata
        array_index += total_parameters_identified
    return all_identifiability_values, number_of_parameters_per_flux


def ident_parameter_name(parameter_id):
    parameter_list = ['V1max', 'K1ac', 'k1cat', 'K1ac',
                      'V2max', 'K2pep',
                      'V3max_option1', 'K3fdp_option1', 'K3pep_option1',
                      'V3max_option2', 'K3fdp_option2', 'K3pep_option2']
    try:
        return [parameter_list[id] for id in parameter_id]
    except TypeError:
        return parameter_list[parameter_id]


def kotte_parameter_name(parameter_id):
    parameter_list = ['K1ac', 'K3fdp', 'L3fdp', 'K3pep',
                      'K2pep', 'vemax', 'Kefdp', 'ne', 'd',
                      'V4max', 'k1cat', 'V3max', 'V2max', 'ac']
    try:
        return [parameter_list[id] for id in parameter_id]
    except TypeError:
        return parameter_list[parameter_id]


def parameter_based_processing(ident_details):
    """parameter-based classification of experimental datasets"""
    number_data, p = ident_details["boolean"].shape
    # p = np.cumsum([0] + parameters_per_flux).tolist()
    fp_list_keys = ['flux{}p{}'.format(f_index + 1, p_index + 1)
                    for f_index, p_limit in enumerate(ident_details["parameters"])
                    for p_index in range(0, p_limit)]
    fp_list_values = [[data_id for data_id in it.compress(range(0, number_data), list(ident_details["boolean"][:, parameter_id]))]
                      for parameter_id in range(0, p)]
    fp_list = dict(zip(fp_list_keys, fp_list_values))
    return fp_list


def data_based_processing(ident_details):
    """dataset dependent classification of parameters"""
    number_data, p = ident_details["boolean"].shape
    data_list_keys = ['dataset{}'.format(perturbation_id + 1) for perturbation_id in range(0, number_data)]
    data_list_values = [[parameter_id for parameter_id in it.compress(range(0, p), list(ident_details["boolean"][data_id, :]))]
                        for data_id in range(0, number_data)]
    data_list = dict(zip(data_list_keys, data_list_values))
    return data_list


def process_ident_data(ident_values, number_data):
    # create signed boolean array for identifiability
    signed_ident_values = np.sign(ident_values)
    ident_fun_val = []
    for id in range(0, number_data):
        ident_fun_val.append(signed_ident_values[id * 12:(id + 1) * 12, -1])
    p_list = [[p_id for p_id, val in enumerate(data_set) if val > 0] for data_set in ident_fun_val]
    p_list_boolean = [[True if parameter_id in list_1 else False for parameter_id in range(0, 12)] for list_1 in p_list]
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


def get_most_useful_dataset(ident_boolean_array):
    """get data identifying most parameters (most useful data set)"""
    # most useful dataset - based on number of parameter identified
    number_data, number_parameter = ident_boolean_array.shape
    identified_parameters = []  # np.zeros((number_data, 1))
    for idata in range(0, number_data):
        ip = [p_id for p_id in it.compress(range(0, number_parameter), ident_boolean_array[idata, :])]
        identified_parameters.append(len(ip))

    # get data identifying most parameters (most useful data set)
    max_useful_data = max(identified_parameters)
    max_data_id = [identified_parameters.index(max(identified_parameters))]
    for idata, number_identified_parameter in enumerate(identified_parameters):
        if number_identified_parameter == max_useful_data and idata not in max_data_id:
            max_data_id.append(idata)
        elif number_identified_parameter > max_useful_data:
            max_useful_data = number_identified_parameter
            max_data_id = [idata]
    max_data = {"maximum": max_useful_data,
                "id": max_data_id,
                "parameters": identified_parameters}
    return max_data


def get_useful_data_info(ident_boolean_array, experiment_details, perturbation_details, data_needed):
    """get info all data sets that identify maximum number of parameters"""
    _, p = ident_boolean_array.shape
    # data_needed = max_data["id"]
    parameters_identified_by_max_data_id = []
    # get parameters identified by data_needed
    data_needed_id = 0
    for data_id in data_needed:
        data_set_name = 'option {}'.format(data_needed_id + 1)
        data_set_ident = [ident_parameter_name(i)
                          for i in it.compress(range(0, p), list(ident_boolean_array[data_id, :]))]
        data_set_ident_id = [i for i in it.compress(range(0, p), list(ident_boolean_array[data_id, :]))]
        # identify all experiments for given data set
        perturbation_id = [int(i) for i in experiment_details["details"][data_id, range(0, 12, 4)].tolist()]
        # identify parameters perturbed by said experiment
        perturbed_parameters = perturbation_details["indices"][perturbation_id, :]
        # get parameter names from parameter indices
        perturbed_parameter_name = [(kotte_parameter_name(int(j[0])), j[1]) for j in perturbed_parameters.tolist()]
        info = {'id': data_id, 'parameters': data_set_ident,
                'experiments': perturbed_parameter_name, 'parameter_ids':data_set_ident_id}
        parameters_identified_by_max_data_id.append((data_set_name, info))
        data_needed_id += 1
    data_list = dict(parameters_identified_by_max_data_id)
    return data_list


def get_additional_data(ident_details, experiment_details, perturbation_details, unidentified_parameter_ids):
    """find data sets that identify parameters in unidentified set"""
    number_data, p = ident_details["boolean"].shape
    # get boolean for unidentified parameters
    unident_boolean = ident_details["boolean"][:, unidentified_parameter_ids]
    # work on this later
    # max_data = get_most_useful_dataset(unident_boolean)
    # if max_data["maximum"] > 1:
    #    selected_data_id = max_data["id"]
    #else:

    # get every single data set that identifies unidentified parameters
    data_ident = [[k for k in it.compress(range(0, number_data), list(unident_boolean[:, unid_p]))]
                  for unid_p in range(0, len(unidentified_parameter_ids))]
    # get other info for data sets in this list
    extra_data_details = []
    for list_elements in data_ident:
        if list_elements:
            data_list = get_useful_data_info(ident_details["boolean"],
                                             experiment_details,
                                             perturbation_details,
                                             list_elements)
            extra_data_details.append(data_list)
        else:
            extra_data_details.append({})



    # collect common data ids in data_ident
    #for j_data_ident in data_ident:
    #for i_data in range(0, number_data):
    #   unidentifiable_parameters = [k for k in it.compress(range(0, p), list(np.logical_not(ident_details["boolean"][i_data, :])))]

    #    identified_parameters = [i for i in it.compress(range(0, p), list(ident_details["boolean"][i_data, :]))]
    #    identifiable_parameter_boolean = [True if un_id in identified_parameters
    #                                      else False
    #                                      for un_id in unidentified_parameter_ids]
    #    if any(identifiable_parameter_boolean):
    #        identifiable_parameter = [j
    #                                  for j in it.compress(unidentified_parameter_ids, identifiable_parameter_boolean)]
    #        extra_data_details.append({'data':i_data, 'parameter':identifiable_parameter})
    return extra_data_details


def calculate_experiment_combos(ident_details, experiment_details, perturbation_details, data_list):
    """choose additional data sets and consequently, experiments (if possible) to identify other parameters not
    identified by chosen data set(s)"""
    # get parameters identified by chosen data set
    number_data, p = ident_details["boolean"].shape
    new_combos = []
    new_combo_number = 0
    for option_id in data_list:
        chosen_data = data_list[option_id]
        data_set_id = chosen_data["id"]
        # get parameters identified
        identified_parameters = [i for i in it.compress(range(0, p), list(ident_details["boolean"][data_set_id, :]))]
        # get unidentified parameters
        nonidentified_parameters = list(set(range(0, p)) - set(identified_parameters))

        # go through all data sets to identify indexes that identify unidentified parameters
        extra_data = get_additional_data(ident_details,
                                         experiment_details,
                                         perturbation_details,
                                         nonidentified_parameters)
        for extra_id, dict_id in enumerate(extra_data):
            for id, options in enumerate(dict_id):
                new_option_name = 'combination {}'.format(new_combo_number + 1)
                new_option_composition = [chosen_data["id"], dict_id[options]["id"]]
                new_option_experiments = [chosen_data["experiments"], dict_id[options]["experiments"]]
                new_option_parameters = list(set(chosen_data["parameters"]) | set(dict_id[options]["parameters"]))
                new_options_parameter_ids = list(
                    set(chosen_data["parameter_ids"]) | set(dict_id[options]["parameter_ids"]))
                combo_info = {'id': new_option_composition,
                              'experiments': new_option_experiments, 'parameters': new_option_parameters,
                              'parameter_ids': new_options_parameter_ids}
                new_combos.append((new_option_name, combo_info))
                new_combo_number += 1
    return new_combos


def process_info(ident_details, experiment_details, perturbation_details, number_fluxes):
    number_data, p = ident_details["boolean"].shape
    # parameter-based classification of experimental datasets
    fp_list = parameter_based_processing(ident_details)

    # dataset dependent classification of parameters
    # data_list = data_based_processing(ident_details)

    # most useful dataset - based on number of parameter identified
    max_data = get_most_useful_dataset(ident_details["boolean"])

    # decide which experiments to perform for each parameter based on above calculations
    # get all info on most identifiable data sets (most useful data set)
    data_list = get_useful_data_info(ident_details["boolean"],
                                     experiment_details,
                                     perturbation_details,
                                     max_data["id"])

    # choose additional data sets and consequently, experiments (if possible) to identify other parameters not
    # identified by chosen data set(s)
    new_combos = calculate_experiment_combos(ident_details, experiment_details, perturbation_details, data_list)

    # most easily identifieable parameter - based on frequency of identification
    identifying_data = [] # np.zeros((number_parameter, 1))
    for iparameter in range(0, p):
        idata = [data_id for data_id in it.compress(range(0, number_data), ident_details["boolean"][:, iparameter])]
        identifying_data.append(len(idata))

    # get parameters identified by most data sets (most identifiable parameter)
    max_identified = max(identifying_data)
    max_parameter_id = [identifying_data.index(max(identifying_data))]
    for iparameter, number_identifying_data in enumerate(identifying_data):
        if number_identifying_data == max_identified and iparameter not in max_parameter_id:
            max_parameter_id.append(iparameter)
        elif number_identifying_data > max_identified:
            max_identified = number_identifying_data
            max_parameter_id = [iparameter]
    max_parameter = {"maximum":max_identified,
                     "id":max_parameter_id,
                     "data":identifying_data}

    # ident_parameter_names = ident_parameter_name(range(0, 12))

    return data_list, new_combos, max_parameter


def establish_kotte_flux_identifiability(all_data, choose):
    """call all identifiability evaluation funcs above and print numerical results"""
    chosen_values = list(all_data["values"][0:choose, :])
    # chosen_details = all_data["details"][0:choose, :]

    ident_fun_list = (flux_1_ident_expression, flux_2_ident_expression, flux_3_ident_expression)
    # number_fluxes = len(ident_fun_list)
    number_data = len(chosen_values)
    ident_values, parameters_per_flux = get_ident_value(ident_fun_list, chosen_values)

    # process identifiability data
    _, plist_boolean = process_ident_data(ident_values, number_data)
    ident_details = {"boolean":plist_boolean, "values":ident_values, "parameters":parameters_per_flux}

    return ident_details


def parameter_change(new_value, old_value):
    return (new_value-old_value)/old_value


def get_changed_parameters(original_parameters, changed_parameters, experiment_index, parameter_index):
    parameter_array = np.zeros((1, 4))
    parameter_array[0, 0] = experiment_index
    parameter_array[0, 1] = parameter_index
    parameter_array[0, 2] = changed_parameters[experiment_index, parameter_index]
    parameter_array[0, 3] = parameter_change(changed_parameters[experiment_index, parameter_index],
                                              original_parameters[parameter_index])
    return parameter_array


def select_experiment_for_dataset(experiment_id):
    return None

def loop_through_experiments(experiments_per_dataset, total_experiments):
    chosen_ids = [0]
    for experiment in experiments_per_dataset:

        pass
    experimental_data = []
    return experimental_data


def arrange_experimental_data(xss, fss,
                              perturbation_details, experiments_per_set, flux_id=np.array([0, 1, 2, 3, 4, 5])):
    """get combinations of datasets using xss, fss and parameters for kotte model.
    see identifiability functions above for order of values in datasets"""
    parameters = perturbation_details["values"]
    experiment_indices = perturbation_details["indices"]
    original = perturbation_details["original"]

    # get all experiment combinations based on experiments_per_set
    number_of_experiments = len(xss)
    data_combinations = []
    for item in it.permutations(range(0, number_of_experiments), experiments_per_set):
        data_combinations.append(item)

    data_combination_boolean = [[True if experiment_id in list_1 else False
                                 for experiment_id in range(0, number_of_experiments)]
                                for list_1 in data_combinations]
    data_combination_boolean = np.array(data_combination_boolean)

    # initialize vectors/arrays and other lists for arrangement
    number_data_sets = len(data_combinations)
    size_of_data_set = 8
    dataset_value_array = np.zeros((number_data_sets, size_of_data_set * experiments_per_set))
    all_parameter_changes = np.zeros((number_data_sets, 4 * experiments_per_set))

    # arrange and collect parameter/ss values
    for data_index, experiment_index in enumerate(data_combinations):
        data = []
        changes = []
        for index in experiment_index:
            data.append(np.hstack((parameters[index][-1],
                                   xss[index],
                                   fss[index][flux_id])))
            changes.append(get_changed_parameters(original_parameters=original,
                                                  changed_parameters=parameters,
                                                  experiment_index=index,
                                                  parameter_index=int(experiment_indices[index, 0])))
        dataset_value_array[data_index, :] = np.hstack(data[:])
        all_parameter_changes[data_index, :] = np.hstack(changes[:])

    dataset_keys = ['values', 'details', 'boolean']
    experimental_data = dict(zip(dataset_keys, [dataset_value_array, all_parameter_changes, data_combination_boolean]))
    return experimental_data
