import numpy as np
from sympy import *
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
    all_flux = np.hstack((flux_1,flux_2,flux_3,flux_4,flux_5,flux_6))

    return all_flux


def kotte_ode(t, y, par_val):

    # K1ac, K3fdp, L3fdp, K3pep, K2pep, vemax, Kefdp, ne, d, V4max, k1cat, V3max, V2max, ac = \
    #     [.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, .1]
    # par_val = np.vstack((K1ac, K3fdp, L3fdp, K3pep, K2pep, vemax, Kefdp, ne, d, V4max, k1cat, V3max, V2max, ac))

    flux = kotte_flux(y,par_val)
    yd_pep = flux[0] - flux[3] - flux[4]
    yd_fdp = flux[3] - flux[2]
    yd_e = flux[1] - flux[5]

    return np.hstack((yd_pep,yd_fdp,yd_e))


def define_sym_variables():
    """define all required symbolic variables for sympy expressions"""
    ac1, ac2, ac3, x11, x12, x13, x21, x22, x23, x31, x32, x33, \
    v31, v32, v33, v11, v12, v13, v21, v22, v23, v41, v42, v43 = \
        symbols('ac1, ac2, ac3, x11, x12, x13, x21, x22, x23, x31, x32, x33,'
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

    # symbolic expression for flux v1 w/o enzyme concentration data
    v1max_numerator = ac1 * v11 * v12 - ac2 * v11 * v12
    v1max_denominator = -(ac2 * v11 - ac1 * v12)
    k1ac_numerator = ac1 * (ac2 * v11 - ac2 * v12)
    k1ac_denominator = -ac2 * v11 + ac1 * v12
    v1max_nr_expr = lambdify([variables], v1max_numerator, "numpy")
    v1max_dr_expr = lambdify([variables], v1max_denominator, "numpy")
    k1ac_nr_expr = lambdify([variables], k1ac_numerator, "numpy")
    k1ac_dr_expr = lambdify([variables], k1ac_denominator, "numpy")
    v1max_expression = lambdify([variables], v1max_numerator/v1max_denominator, "numpy")
    k1ac_expression = lambdify([variables], k1ac_numerator/k1ac_denominator, "numpy")
    v1max_no_enzyme_numerator_value = v1max_nr_expr(experimental_data)
    v1max_no_enzyme_denominator_value = v1max_dr_expr(experimental_data)
    v1max_no_enzyme_value = v1max_expression(experimental_data)
    k1ac_no_enzyme_numerator_value = k1ac_nr_expr(experimental_data)
    k1ac_no_enzyme_denominator_value = k1ac_dr_expr(experimental_data)
    k1ac_no_enzyme_value = k1ac_expression(experimental_data)

    # symbolic expression for flux v1 w/ enzyme concentration data
    k1cat_nr = - ac1 * v11 * v12 + ac2 * v11 * v12
    k1cat_dr = -(ac1 * v12 * x31 - ac2 * v11 * x32)
    k1cat_nr_expr = lambdify([variables], k1cat_nr, "numpy")
    k1cat_dr_expr = lambdify([variables], k1cat_dr, "numpy")
    k1cat_expression = lambdify([variables], k1cat_nr / k1cat_dr, "numpy")
    k1ac_nr = ac1 * (-ac2 * v12 * x31 + ac2 * v11 * x32)
    k1ac_dr = ac1 * v12 * x31 - ac2 * v11 * x32
    k1ac_nr_expr = lambdify([variables], k1ac_nr, "numpy")
    k1ac_dr_expr = lambdify([variables], k1ac_dr, "numpy")
    k1ac_expression = lambdify([variables], k1ac_nr/k1ac_dr, "numpy")
    k1cat_enzyme_numerator_value = k1cat_nr_expr(experimental_data)
    k1cat_enzyme_denominator_value = k1cat_dr_expr(experimental_data)
    k1cat_enzyme_value = k1cat_expression(experimental_data)
    k1ac_enzyme_numerator_value = k1ac_nr_expr(experimental_data)
    k1ac_enzyme_denominator_value = k1ac_dr_expr(experimental_data)
    k1ac_enzyme_value = k1ac_expression(experimental_data)

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
    v2max_nr_expr = lambdify([variables], v2max_numerator, "numpy")
    v2max_dr_expr = lambdify([variables], v2max_denominator, "numpy")
    v2max_expr = lambdify([variables], v2max_numerator/v2max_denominator, "numpy")
    k2pep_nr_expr = lambdify([variables], k2pep_numerator, "numpy")
    k2pep_dr_expr = lambdify([variables], k2pep_denominator, "numpy")
    k2pep_expr = lambdify([variables], k2pep_numerator/k2pep_denominator, "numpy")
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
    v3max_numerator = -v31*v32*v33*(x11*x12*x21 + x11*x13*x21 + x11*x12*x22 -
                                    x12*x13*x22 - x11*x13*x23 + x12*x13*x23) - \
                    (v31*v32*v33*x12*x21*(v33*x21*(-v31*x11*x12*x22 + v32*x11*x12*x22 + v32*x11*x13*x22 -
                                                   v31*x12*x13*x22 - v32*x11*x12*x23 - v32*x11*x13*x23)  +
                    v31*x23*(v32*x11*x13*x21 + v32*x12*x13*x21 + v33*x11*x12*x22 -
                             v32*x11*x13*x22 - v32*x12*x13*x22 + v33*x12*x13*x22) -
                    sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 - v32*v33*x11*x13*x21*x22 +
                          v31*v33*x12*x13*x21*x22 + v32*v33*x11*x12*x21*x23 - v31*v32*x11*x13*x21*x23 +
                          v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x21*x23 - v31*v33*x11*x12*x22*x23 +
                          v31*v32*x11*x13*x22*x23 + v31*v32*x12*x13*x22*x23 - v31*v33*x12*x13*x22*x23)**2 -
                         4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                            v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                         (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 - v31*v32*x11*x12*x13*x21*x23 +
                          v32*v33*x11*x12*x13*x21*x23 + v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))))/ \
                    (2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 -
                        v31*v32*x13*x21*x23 - v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)) + \
                    (v31*v32*v33*x13*x21*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                          v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                          v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                          v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                          v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                          v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 -
                                          sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
                                                v32*v33*x11*x13*x21*x22 + v31*v33*x12*x13*x21*x22 +
                                                v32*v33*x11*x12*x21*x23 - v31*v32*x11*x13*x21*x23 +
                                                v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x21*x23 -
                                                v31*v33*x11*x12*x22*x23 +  v31*v32*x11*x13*x22*x23 +
                                                v31*v32*x12*x13*x22*x23 - v31*v33*x12*x13*x22*x23)**2 -
                                               4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 +
                                                  v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                                                  v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                                               (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 -
                                                v31*v32*x11*x12*x13*x21*x23 + v32*v33*x11*x12*x13*x21*x23 +
                                                v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))))/ \
                    (2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 -
                        v31*v32*x13*x21*x23 - v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)) + \
                    (v31*v32*v33*x11*x22*(-v31*v33*x11*x12*x21*x22 - v31*v33*x12*x13*x21*x22 +
                                          v31*v33*x12*x13*x22*x23 + v31*v33*x11*x12*x22*x23 +
                                          v32*v33*x11*x12*x21*x22 + v32*v33*x11*x13*x21*x22 -
                                          v32*v33*x11*x12*x21*x23 - v32*v33*x11*x13*x21*x23 +
                                          v31*v32*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 -
                                          v31*v32*x11*x13*x22*x23 - v31*v32*x12*x13*x22*x23  -
                                          sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
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
                    (2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 -
                        v31*v32*x13*x21*x23 - v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)) - \
                    (v31*v32*v33*x13*x22*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x21*x22*(x12 + x13) -
                                          v31*v33*x12*x13*x21*x22 - v32*v33*x11*x21*x23*(x12 + x13) +
                                          v31*v32*x13*x21*x23*(x11 + x12) + v31*v33*x12*x22*x23*(x11 + x13) -
                                          v31*v32*x13*x22*x23*(x11 + x12) -
                                          sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
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
                    (2*(v31*v33*x12*x22*(x21 - x23) - v32*v33*x11*x21*(x22 - x23) - v31*v32*x13*x23*(x21 + x22))) - \
                    (v31*v32*v33*x11*x23 (-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                          v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                          v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                          v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                          v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                          v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 -
                                          sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11**x21*x22(x12 - x13) +
                                                v31*v33*x12*x13*x21*x22 + v32*v33*x11*x21*x23(x12 + x13) -
                                                v31*v32*x13*x21*x23*(x11 + x12) - v31*v33*x12*x22*x23*(x11 + x13) +
                                                v31*v32*x13*x22*x23*(x11 + x12))**2 -
                                               4*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 +
                                                  v32*v33*x11*x21*x23 - v31*v32*x13*x21*x23 -
                                                  v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)*
                                               (v31*v33*x11*x12*x13*x21*x22 - v32*v33*x11*x12*x13*x21*x22 -
                                                v31*v32*x11*x12*x13*x21*x23 + v32*v33*x11*x12*x13*x21*x23 +
                                                v31*v32*x11*x12*x13*x22*x23 - v31*v33*x11*x12*x13*x22*x23))))/\
                    (2*(-v32*v33*x11*x21*x22 + v31*v33*x12*x21*x22 + v32*v33*x11*x21*x23 -
                        v31*v32*x13*x21*x23 - v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)) + \
                      (v31*v32*v33*x12*x23*(-v31*v33*x12*x22*(x11*x21 + x13*x21 - x11*x23 - x13*x23) +
                                            v32*v33*x11*x21*(x12*x22 + x13*x22 - x12*x23 - x13*x23) +
                                            v31*v32*x13*x23*(x11*x21 + x12*x21 - x11*x22 - x12*x22) -
                                            sqrt((v31*v33*x12*x22*(x11*x21 + x13*x21 - x11*x23 - x13*x23) -
                                                  v32*v33*x11*x21*(x22*x12 + x13*x22 - x12*x23 - x13*x23) -
                                                  v31*v32*x13*x23*(x11*x21 - x12*x21 - x11*x22 - x12*x22))**2 -
                                                 4*(v31*v33*x12*x22*(x21 - x23) - v31*v32*x13*x23*(x21 - x22) -
                                                    v32*v33*x11*x21*(x22 - x23))*(v31*v33*x11*x12*x13*x22*(x21 - x23) -
                                                                                  v32*v33*x11*x12*x13*x21*(x22 - x23) -
                                                                                  v31*v32*x11*x12*x13*x23(x21 - x22)))))/\
                    (2*(-v32*v33*x11*x21*(x22 - x23) + v31*v33*x12*x22*(x21 - x23) - v31*v32*x13*x23*(x21 - x22)))
    v3max_nr_expr = lambdify([variables], v3max_numerator, "numpy")
    v3max_nr_value = v3max_nr_expr(experimental_data)
    v3max_denominator = -v32 * v33 * x11 * x12 * x21 + v32 * v33 * x11 * x13 * x21 + v31 * v33 * x11 * x12 * x22 - \
                  v31 * v33 * x12 * x13 * x22 - v31 * v32 * x11 * x13 * x23 + v31 * v32 * x12 * x13 * x23
    v3max_dr_expr = lambdify([variables], v3max_denominator, "numpy")
    v3max_dr_value = v3max_dr_expr(experimental_data)
    v3max_expr = lambdify([variables], v3max_numerator/v3max_denominator, "numpy")
    v3max_value = v3max_expr(experimental_data)
    # K3fdp
    v3max_numerator_2 = -v31*v32*v33*(x11*x12*x21 - x11*x13*x21 - x11*x12*x22 + x12*x13*x22 + x11*x13*x23 - x12*x13*x23) - \
                    (v31*v32*v33*x12*x21*(-v31*v33*x12*x22*(x11*x21  + x13*x21 - x11*x23 - x13*x23) +
                                            v32*v33*x11*x21*(x12*x22 + x13*x22 - x12*x23 - x13*x23) +
                                          v31*v32*x13*x23*(x11*x21 + x12*x21 - x11*x22 - x12*x22)  +
                                          sqrt((v31*v33*x12*x22*(x11*x21 + x13*x21 - x11*x23 - x13*x23) -
                                                v32*v33*x11*x21*(x12*x22 + x13*x22 - x12*x23 - x13*x23) -
                                                v31*v32*x13*x23*(x12*x21 - x11*x22 - x12*x22 + x11*x21))**2 -
                                               4*(-v32*v33*x11*x21*(x22 - x23) - v31*v32*x13*x23*(x21 - x22) -
                                                  v31*v33*x12*x22*(x23 + x21))*(v31*v33*x11*x12*x13*x22*(x21 - x23) -
                                                                                v32*v33*x11*x12*x13*x21*(x22 - x23) -
                                                                                v31*v32*x11*x12*x13*x23*(x21 - x22)))))/\
                    (2*(-v32*v33*x11*x21*(x22 - x23) + v31*v33*x12*x22*(x21 - x23) - v31*v32*x13*x23*(x21 - x22))) + \
                        (v31*v32*v33*x13*x21*(v31*v33*x12*x22*(x11*x23 + x13*x23 - x11*x21 - x13*x21) +
                                          v32*v33*x11*x21*(x12*x22 + x13*x22 - x12*x23 - x13*x23) +
                                          v31*v32*x13*x23*(x11*x21 + x12*x21 - x11*x22 - x12*x22) +
                                          sqrt((v31*v33*(x11*x12*x21*x22 + x12*x13*x21*x22 -
                                                         x12*x13*x22*x23 - x11*x12*x22*x23) -
                                                v32*v33*(x11*x12*x21*x22 + x11*x13*x21*x22  -
                                                         x11*x12*x21*x23 - x11*x13*x21*x23) -
                                                v31*v32*(x11*x13*x21*x23 + x12*x13*x21*x23 -
                                                         x11*x13*x22*x23 - x12*x13*x22*x23))**2 -
                                               4*(v32*v33*(x11*x21*x23 - x11*x21*x22) +
                                                  v31*v33*(x12*x21*x22 - x12*x22*x23) +
                                                  v31*v32*(x13*x22*x23 - x13*x21*x23))*
                                               (v31*v33*(x11*x12*x13*x21*x22 - x11*x12*x13*x22*x23) -
                                                v32*v33*(x11*x12*x13*x21*x22 - x11*x12*x13*x21*x23) +
                                                v31*v32*(x11*x12*x13*x22*x23 - x11*x12*x13*x21*x23)))))/\
                        (2*(-v32*v33*(x11*x21*x22 - x11*x21*x23) + v31*v32*(x13*x22*x23 - x13*x21*x23) +
                            v31*v33*(x12*x21*x22 - x12*x22*x23))) +\
                        (v31*v32*v33*x11*x22*(-v32*v33*(x11*x12*x21*x23 + x11*x13*x21*x23 -
                                                        x11*x12*x21*x22 - x11*x13*x21*x22) +
                                              v31*v32*(x11*x13*x21*x23 + x12*x13*x21*x23 -
                                                       x11*x13*x22*x23 - x12*x13*x22*x23) +
                                              v31*v33*(x11*x12*x22*x23 + x12*x13*x22*x23 -
                                                       x11*x12*x21*x22 - x12*x13*x21*x22) +
                                              sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
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
                            v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)) -\
                        (v31*v32*v33*x13*x22*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                              v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                              v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                              v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                              v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                              v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 +
                                              sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
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
                        (v31*v32*v33*x11*x23*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                              v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                              v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                              v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                              v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                              v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 +
                                              sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
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
                            v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23)) + \
                        (v31*v32*v33*x12*x23*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                              v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                              v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                              v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                              v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                              v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 +
                                              sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
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
                            v31*v33*x12*x22*x23 + v31*v32*x13*x22*x23))
    k3fdp_numerator = -v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 + v31*v32*x11*x13*x21*x23 - \
                      v32*v33*x11*x13*x21*x23 - v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 + \
                      (v32*v33*x11*x21*x22*(-v31*v33*x11*x12*x21*x22 + v32*v33*x11*x12*x21*x22 +
                                            v32*v33*x11*x13*x21*x22 - v31*v33*x12*x13*x21*x22 -
                                            v32*v33*x11*x12*x21*x23 + v31*v32*x11*x13*x21*x23 -
                                            v32*v33*x11*x13*x21*x23 + v31*v32*x12*x13*x21*x23 +
                                            v31*v33*x11*x12*x22*x23 - v31*v32*x11*x13*x22*x23 -
                                            v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 -
                                            sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
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
                                            v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 -
                                            sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
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
                                            v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 -
                                            sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
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
                                            v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 -
                                            sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
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
                                            v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 -
                                            sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
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
                                            v31*v32*x12*x13*x22*x23 + v31*v33*x12*x13*x22*x23 -
                                            sqrt((v31*v33*x11*x12*x21*x22 - v32*v33*x11*x12*x21*x22 -
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
    k3fdp_denominator = -v32 * v33 * x11 * x12 * x21 + v32 * v33 * x11 * x13 * x21 + v31 * v33 * x11 * x12 * x22 - \
                  v31 * v33 * x12 * x13 * x22 - v31 * v32 * x11 * x13 * x23 + v31 * v32 * x12 * x13 * x23
    # K3pep
    k3pep_denominator = 2 * (-v32 * v33 * x11 * x21 * x22 + v31 * v33 * x12 * x21 * x22 + v32 * v33 * x11 * x21 * x23 -
                       v31 * v32 * x13 * x21 * x23 - v31 * v33 * x12 * x22 * x23 + v31 * v32 * x13 * x22 * x23)

    # K3fdp_sol_2
    # K3pep_sol_2
    # V3max_sol_2 = -v32*v33*x11*x12*x21 + v32*v33*x11*x13*x21 + v31*v33*x11*x12*x22 - \
    #               v31*v33*x12*x13*x22 - v31*v32*x11*x13*x23 + v31*v32*x12*x13*x23


    k3fdp_dr_expr = lambdify([variables], k3fdp_denominator, "numpy")
    k3pep_dr_expr = lambdify([variables], k3pep_denominator, "numpy")

    k3fdp_dr_value = k3fdp_dr_expr(experimental_data)
    k3pep_dr_value = k3pep_dr_expr(experimental_data)

    return [v3max_nr_value, v3max_dr_value, v3max_value, k3fdp_dr_value, k3pep_dr_value]
