import casadi as casadi
import numpy as np
from kotte_model import ident_parameter_name
from kotte_model import kotte_true_parameter_values


def v3_ck_numerical_problem(chosen_data):
    """fun to calculate objectives and constraints for numerical identification of v3
    using convenience kinetics"""

    number_data = len(chosen_data)
    # symbolic variables for use with casadi
    v3max = casadi.SX.sym("v3max")
    k3fdp = casadi.SX.sym("k3fdp")
    k3pep = casadi.SX.sym("k3pep")
    error = casadi.SX.sym("error", number_data)

    all_constraints = []
    for i_data_id, i_data in enumerate(chosen_data):
        _, pep, fdp, _, _, _, v3, _, _, _ = list(i_data)

        # flux equation for each experiment i_data in chosen_data
        fdp_sat = fdp / k3fdp
        pep_sat = pep / k3pep
        nr_3 = v3max * fdp_sat
        dr_3 = 1 + fdp_sat
        regulation_activate = 1 / (1 + 1 / pep_sat)
        flux_3 = regulation_activate * nr_3 / dr_3
        # formulate constraint for each experimental data set
        all_constraints.append(flux_3 - v3 - error[i_data_id])

    # create casadi array of constraints
    constraints = casadi.vertcat(*all_constraints)

    # create complete casadi objective fun
    objective = error[0]**2
    for i_error_id in range(1, number_data):
        objective += error[i_error_id]**2

    # get complete set of all optimization variables vertcat([parameters, error])
    x = casadi.vertcat(*[v3max, k3fdp, k3pep, error])

    nlp = {'x': x, 'f': objective, 'g': constraints}

    return nlp


def v3_numerical_problem(chosen_data):
    """fun to calculate objectives and constraints for numerical identification of v3
        using MWC kinetics"""

    number_data = len(chosen_data)
    # symbolic variables for use with casadi
    v3max = casadi.SX.sym("v3max")
    k3fdp = casadi.SX.sym("k3fdp")
    k3pep = casadi.SX.sym("k3pep")
    l3fdp = casadi.SX.sym("l3fdp")
    error = casadi.SX.sym("error", number_data)

    # fdp_sat = 1 + y[1] / K3fdp
    # pep_sat = 1 + y[0] / K3pep
    # flux_3 = V3max * (fdp_sat - 1) * (fdp_sat ** 3) / (fdp_sat ** 4 + L3fdp * (pep_sat ** (-4)))

    all_constraints = []
    for i_data_id, i_data in enumerate(chosen_data):
        _, pep, fdp, _, _, _, v3, _, _, _ = list(i_data)

        # flux equation for each experiment i_data in chosen_data
        fdp_sat = 1 + fdp / k3fdp
        pep_sat = 1 + pep / k3pep
        nr_3 = v3max * (fdp_sat - 1) * (fdp_sat**3)
        dr_3 = (fdp_sat ** 4 + l3fdp * 1e6 * pep_sat ** (-4))
        flux_3 = nr_3 / dr_3
        # formulate constraint for each experimental data set
        all_constraints.append(flux_3 - v3 - error[i_data_id])

    # create casadi array of constraints
    constraints = casadi.vertcat(*all_constraints)

    # create complete casadi objective fun
    objective = error[0] ** 2
    for i_error_id in range(1, number_data):
        objective += error[i_error_id] ** 2

    # get complete set of all optimization variables vertcat([parameters, error])
    x = casadi.vertcat(*[v3max, k3fdp, k3pep, l3fdp, error])

    nlp = {'x': x, 'f': objective, 'g': constraints}

    return nlp


def solve_numerical_nlp(chosen_fun, chosen_data, opt_problem_details, optim_options={}):
    """solve nlp for determining parameter values numerically"""
    # formulate nlp
    nlp = chosen_fun(chosen_data)

    # NLP solver options
    if optim_options:
        solver = optim_options["solver"]
        opts = optim_options["opts"]
    else:
        solver = "ipopt"
        opts = {"ipopt.tol": 1e-12}
        # solver = "sqpmethod"
        # opts = {"qpsol": "qpoases"}

    # Allocate an NLP solver and buffer
    solver = casadi.nlpsol("solver", solver, nlp, opts)

    # Solve the problem
    res = solver(**opt_problem_details)

    return res


def parse_opt_result(opt_sol):
    """convert all optimization results from casadi.DM to numpy arrays"""
    new_opt_sol = {"f": opt_sol["f"].full(),
                   "xopt": opt_sol["x"].full(),
                   "lam_x": opt_sol["lam_x"].full(),
                   "lam_g": opt_sol["lam_g"].full()}
    # Print the optimal cost
    print("optimal cost: ", float(new_opt_sol["f"]))

    # Print the optimal solution
    print("optimal solution: ", new_opt_sol["xopt"])
    print("dual solution (x) = ", new_opt_sol["lam_x"])
    print("dual solution (g) = ", new_opt_sol["lam_g"])
    return new_opt_sol


def opt_result_for_plots(all_data_set_opt_sol):
    """parse optimization results from all data sets to generate data for box or other(?) plot"""
    # number_data = len(all_data_set_opt_sol)

    all_obj_fun = [i_data_sol["f"] for i_data_sol in all_data_set_opt_sol]
    all_xopt = [i_data_sol["xopt"] for i_data_sol in all_data_set_opt_sol]
    all_lam_x = [i_data_sol["lam_x"] for i_data_sol in all_data_set_opt_sol]
    all_lam_g = [i_data_sol["lam_g"] for i_data_sol in all_data_set_opt_sol]

    all_solution_details = {"f": all_obj_fun,
                            "x": all_xopt,
                            "lam_x": all_lam_x,
                            "lam_g": all_lam_g}
    return all_solution_details


def identify_all_data_sets(experimental_data, chosen_fun, x0, optim_options={}):
    if chosen_fun == 0:
        ident_fun = v3_ck_numerical_problem
        # Initial condition
        arg = {"lbx": 6 * [0],
               "ubx": [2, 1, 1, .1, .1, .1],
               "lbg": 3 * [0],
               "ubg": 3 * [0]}
    elif chosen_fun == 1:
        ident_fun = v3_numerical_problem
        arg = {"lbx": 8 * [0],
               "ubx": [2, 1, 1, 5, .1, .1, .1, .1],
               "lbg": 4 * [0],
               "ubg": 4 * [0]}
    else:
        ident_fun = []
        arg = {}
    arg["x0"] = x0

    number_data_sets = len(experimental_data)
    all_data_solutions = []
    for i_data_id, i_data in enumerate(experimental_data):
        print("Performing Identifiability Analysis on Data set {} of {}".format(i_data_id+1, number_data_sets))
        if ident_fun:
            sol = solve_numerical_nlp(chosen_fun=ident_fun, chosen_data=i_data, opt_problem_details=arg,
                                      optim_options=optim_options)
        else:
            sol = []
        if sol:
            sol = parse_opt_result(sol)
        all_data_solutions.append(sol)
        print("Identifiability analysis on data set {} of {} complete".format(i_data_id+1, number_data_sets))

    opt_solution = opt_result_for_plots(all_data_solutions)
    return opt_solution


def process_opt_solution(opt_solution, number_of_parameters, flux_id, flux_choice):
    """parse optimization solution based on each parameter from each data set"""
    # optimal values
    x_opt = opt_solution["x"]
    # obj fun value
    obj_f = opt_solution["f"]

    flux_name = ['flux{}'.format(flux_id)] * number_of_parameters
    flux_choice = flux_choice * 3
    all_parameter_name = ident_parameter_name(parameter_id=range(0, number_of_parameters),
                                              flux_name=flux_name, flux_choice_id=flux_choice)
    all_parameter_true_value = kotte_true_parameter_values(flux_based=1, flux_name=flux_name,
                                                           flux_choice_id=flux_choice, parameter_id=all_parameter_name)

    all_parameter_values = []
    for i_parameter in range(number_of_parameters):
        i_parameter_sol = [i_data_set_sol[i_parameter] for i_data_set_sol in x_opt]
        all_parameter_values.append(np.array(i_parameter_sol).flatten())
    all_parameter_info = {"values": all_parameter_values,
                          "names": all_parameter_name,
                          "flux": flux_name,
                          "flux choice": flux_choice,
                          "true value": all_parameter_true_value}
    return all_parameter_info
