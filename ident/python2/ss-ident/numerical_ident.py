import casadi as casadi
import numpy as np
import itertools as it
import pandas as pd
from numpy.random import RandomState
from names_strings import ident_parameter_name
from names_strings import true_parameter_values
from collections import defaultdict
from process_ident_data import write_ident_info_file


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
        _, pep, fdp, _, _, _, v3, _, = list(i_data)

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


def v3_ck_numerical_problem_l1_obj(chosen_data):
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
        _, pep, fdp, _, _, _, v3, _, = list(i_data)

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
    objective = error[0]
    for i_error_id in range(1, number_data):
        objective += error[i_error_id]

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
        opts = {"ipopt.tol": 1e-16}
        # solver = "sqpmethod"
        # opts = {"qpsol": "qpoases"}

    # Allocate an NLP solver and buffer
    solver = casadi.nlpsol("solver", solver, nlp, opts)

    # Solve the problem
    res = solver(**opt_problem_details)

    return res


def parse_opt_result(opt_sol, data_id):
    """convert all optimization results from casadi.DM to numpy arrays"""
    new_opt_sol = {"f": opt_sol["f"].full(),
                   "xopt": opt_sol["x"].full(),
                   "lam_x": opt_sol["lam_x"].full(),
                   "lam_g": opt_sol["lam_g"].full(),
                   "data_set_id": data_id}
    # Print the optimal cost
    print("optimal cost: ", float(new_opt_sol["f"]))

    # Print the optimal solution
    print("optimal solution: ", new_opt_sol["xopt"])
    print("dual solution (x) = ", new_opt_sol["lam_x"])
    print("dual solution (g) = ", new_opt_sol["lam_g"])

    # set optimal
    optimal_bool = [True if abs(i_var) < 1e-14 else False for i_var in new_opt_sol["lam_g"]]
    new_opt_sol["optimal"] = all(optimal_bool)

    return new_opt_sol


def opt_result_for_plots(all_data_set_opt_sol):
    """parse optimization results from all data sets to generate data for box or other(?) plot"""
    # number_data = len(all_data_set_opt_sol)

    all_obj_fun = [i_data_sol["f"] for i_data_sol in all_data_set_opt_sol]
    all_xopt = [i_data_sol["xopt"] for i_data_sol in all_data_set_opt_sol]
    all_lam_x = [i_data_sol["lam_x"] for i_data_sol in all_data_set_opt_sol]
    all_lam_g = [i_data_sol["lam_g"] for i_data_sol in all_data_set_opt_sol]
    all_data_id = [i_data_sol["data_id"] for i_data_sol in all_data_set_opt_sol]

    all_solution_details = {"f": all_obj_fun,
                            "x": all_xopt,
                            "lam_x": all_lam_x,
                            "lam_g": all_lam_g,
                            "data_id": all_data_id}
    return all_solution_details


def compile_opt_result(all_opt_sol):
    """compile optimizaition results to create dict suitable for data frame creation"""
    variable_names = all_opt_sol[0]["variable_name"]
    number_variables = len(variable_names)
    var_id = range(0, number_variables)
    error_names = all_opt_sol[0]["error_name"]
    err_id = [e_id + number_variables for e_id in range(0, len(error_names))]

    all_var_opt = [i_x_value[0] for i_data_sol in all_opt_sol for i_var_id, i_x_value in enumerate(i_data_sol["xopt"])
                   if i_var_id in var_id]
    all_error_opt = [i_x_value[0] for i_data_sol in all_opt_sol for i_var_id, i_x_value in enumerate(i_data_sol["xopt"])
                     if i_var_id in err_id]
    all_var_name = [i_var_name for _ in all_opt_sol for i_var_name in variable_names]
    all_err_name = [i_err_name for _ in all_opt_sol for i_err_name in error_names]
    all_lam_x_names = ["all_lam_x_{}".format(i_var) for i_var in range(0, number_variables + len(error_names))]
    all_var_lam_name = [i_x_value for _ in all_opt_sol for i_var_id, i_x_value in enumerate(all_lam_x_names)
                        if i_var_id in var_id]
    all_err_lam_name = [i_x_value for _ in all_opt_sol for i_var_id, i_x_value in enumerate(all_lam_x_names)
                        if i_var_id in err_id]
    all_var_lam_x = [i_x_value[0] for i_data_sol in all_opt_sol for i_var_id, i_x_value in
                     enumerate(i_data_sol["lam_x"]) if i_var_id in var_id]
    all_err_lam_x = [i_x_value[0] for i_data_sol in all_opt_sol for i_var_id, i_x_value in
                     enumerate(i_data_sol["lam_x"]) if i_var_id in err_id]
    all_cons_lam_name = [i_x_value for i_data_sol in all_opt_sol for i_x_value in i_data_sol["constraint_name"]]
    all_cons_lam_x = [i_x_value[0] for i_data_sol in all_opt_sol for i_x_value in i_data_sol["lam_g"]]
    all_f_opt = [i_obj_value for i_data_sol in all_opt_sol for i_obj_value in
                 [float(i_data_sol["f"][0])] * number_variables]
    all_data_set_id = [i_x_data_set for i_data_sol in all_opt_sol for i_x_data_set in
                       [i_data_sol["data_set_id"]] * number_variables]

    # parameter identifiability based on their numerical values (> 0) use all_var_opt
    all_p_ident = [True if i_var_opt > 1e-4 else False for i_var_opt in all_var_opt]
    # all_lam_opt = [True if abs(i_var_lam) < 1e-16 else False for i_var_lam in all_var_lam_x]
    all_opt = [j_opt for i_sol_opt in all_opt_sol for j_opt in [i_sol_opt["optimal"]] * number_variables]

    all_sol_dict = {"data_set_id": all_data_set_id, "parameter_value": all_var_opt, "parameter_name": all_var_name,
                    "error_value": all_error_opt, "error_name": all_err_name, "variable_dual_name": all_var_lam_name,
                    "error_dual_name": all_err_lam_name, "variable_dual_value": all_var_lam_x,
                    "error_dual_value": all_err_lam_x, "objective": all_f_opt, "identifiability": all_p_ident,
                    "optimality": all_opt, "constraint_dual_name": all_cons_lam_name,
                    "constraint_dual_value": all_cons_lam_x}

    return all_sol_dict


def identify_all_data_sets(experimental_data, chosen_fun, x0, prob, optim_options={}):
    """run numerical identifying algorithm for chosen flux (choose constraints, objective and bounds)
    for a given set of experimental data and initial condition"""
    if chosen_fun == 0:
        ident_fun = v3_ck_numerical_problem
        # set problem constraint and variable bounds
        arg = prob
        variable_name = ["V3max", "K3fdp", "K3pep"]
        error_names = ["eps_1", "eps_2", "eps_3"]
    elif chosen_fun == 1:
        ident_fun = v3_ck_numerical_problem_l1_obj
        # set problem constraint and variable bounds
        arg = prob
        variable_name = ["V3max", "K3fdp", "K3pep"]
        error_names = ["eps_1", "eps_2", "eps_3"]
        cons_name = ["cons_1", "cons_2", "cons_3"]
    elif chosen_fun == 2:
        ident_fun = v3_numerical_problem
        arg = {"lbx": 8 * [0],
               "ubx": [2, 1, 1, 5, .1, .1, .1, .1],
               "lbg": 4 * [0],
               "ubg": 4 * [0]}
        variable_name = ["V3max", "K3fdp", "K3pep", "L3fdp"]
        error_names = ["eps_1", "eps_2", "eps_3", "eps_4"]
        cons_name = ["cons_1", "cons_2", "cons_3", "cons_4"]
    else:
        ident_fun = []
        arg = {}
        variable_name = []
        error_names = []
    # Initial condition
    arg["x0"] = x0

    number_data_sets = len(experimental_data["value"])
    all_data_solutions = []
    for i_data_id, (i_data_value, i_data_name) in enumerate(zip(experimental_data["value"],
                                                                experimental_data["data_set_id"])):
        print("Performing Identifiability Analysis on Data set {} of {}".format(i_data_id+1, number_data_sets))
        if ident_fun:
            sol = solve_numerical_nlp(chosen_fun=ident_fun, chosen_data=i_data_value, opt_problem_details=arg,
                                      optim_options=optim_options)
        else:
            sol = []
        if sol:
            sol = parse_opt_result(sol, i_data_name)
        sol.update({"variable_name": variable_name, "error_name": error_names, "constraint_name": cons_name})
        all_data_solutions.append(sol)
        print("Identifiability analysis on data set {} of {} complete".format(i_data_id+1, number_data_sets))

    opt_solution = compile_opt_result(all_data_solutions)
    # opt_solution = opt_result_for_plots(all_data_solutions)
    return opt_solution


def process_opt_solution(opt_sol, exp_df, opt_solution, number_of_parameters, flux_id, flux_choice):
    """parse optimization solution based on each parameter from each data set"""

    idx = pd.IndexSlice

    # lexicographic ordering of opt_sol df indices
    opt_sol.sort_index(level='sample_name', inplace=True)
    opt_sol.sort_index(level='data_set_id', inplace=True)

    # get all parameter names
    parameter_names = opt_sol["parameter_name"].unique().tolist()

    # get all sample/data set ids (df indices) identifying each parameter
    # sample_data_set_ids = np.unique(opt_sol.index.values).tolist()

    # total number of data sets used
    total_data_sets = [len(opt_sol.index.levels[1])] * len(parameter_names)

    # bring all details together in number_parameter len lists
    relevant_df = opt_sol.loc[:, ['parameter_name', 'parameter_value', 'identifiability', 'optimality']]
    all_parameter_values = []
    all_parameter_names = []
    all_parameter_data_set_ids = []
    for i_parameter, i_parameter_name in enumerate(parameter_names):
        ident_df = relevant_df[(relevant_df["parameter_name"] == i_parameter_name) & (relevant_df["optimality"]) &
                               (relevant_df["identifiability"])]
        # ident_df = optimal_df[(optimal_df["parameter_name"] == i_parameter_name) & (optimal_df["identifiability"])]
        i_p_value = ident_df["parameter_value"].values
        all_parameter_values.append([j_p_value for j_p_value in i_p_value])
        all_parameter_names.append(i_parameter_name)
        # get all sample/data set ids (df indices) identifying each parameter
        all_parameter_data_set_ids.append(ident_df.index.values.tolist())

    # get flux names
    # all_flux_names = ident_df["flux_name"].unique().tolist() * len(all_parameter_info["names"])
    # all_parameter_info.update({"flux_name": all_flux_names})

    all_parameter_info = {"names": all_parameter_names,
                          "values": all_parameter_values,
                          "sample_data_set_id": all_parameter_data_set_ids,
                          "total_data_sets": total_data_sets,
                          "flux_name": ['flux'] * len(parameter_names)}

    # lexographic ordering of exp df indices
    # exp_df.sort_index(level='sample_name', inplace=True)
    # exp_df.sort_index(level='data_set_id', inplace=True)
    # exp_df.sort_index(level='experiment_id', inplace=True)
    #
    #
    # # optimal values
    # x_opt = opt_solution["x"]
    # # obj fun value
    # obj_f = opt_solution["f"]
    #
    # flux_name = ['flux{}'.format(flux_id)] * number_of_parameters
    # flux_choice = flux_choice * 3
    # all_parameter_name = ident_parameter_name(parameter_id=range(0, number_of_parameters),
    #                                           flux_name=flux_name, flux_choice_id=flux_choice)
    # all_parameter_true_value = true_parameter_values(flux_based=1, flux_name=flux_name,
    #                                                        flux_choice_id=flux_choice, parameter_id=all_parameter_name)
    #
    # all_parameter_values = []
    # for i_parameter in range(number_of_parameters):
    #     i_parameter_sol = [i_data_set_sol[i_parameter] for i_data_set_sol in x_opt]
    #     all_parameter_values.append(np.array(i_parameter_sol).flatten())
    # all_parameter_info = {"values": all_parameter_values,
    #                       "names": all_parameter_name,
    #                       "flux": flux_name,
    #                       "flux choice": flux_choice,
    #                       "true value": all_parameter_true_value,
    #                       "data_id": opt_solution["data_id"]}
    return all_parameter_info


def solve_multiple_initial_conditions(all_initial_conditions, experimental_data, chosen_fun, prob, optim_options,
                                      number_of_parameters, flux_id, flux_choice, exp_df, file_name=()):
    """solve numerical nlp ident for multiple parameter initial conditions"""
    number_initial_conditions = len(all_initial_conditions)
    all_x0_all_parameter_opt_info = []
    all_initial_condition_sol = defaultdict(list)
    empty_dict = {}
    for j_id, j_initial_condition in enumerate(all_initial_conditions):
        print("Working on Initial Conditions {} of {}....".format(j_id+1, number_initial_conditions))
        initial_condition_id = "sample_{}".format(j_id)
        # code to run opt for each initial value in each element of list
        opt_solution = identify_all_data_sets(experimental_data, chosen_fun=chosen_fun,
                                              x0=j_initial_condition, prob=prob, optim_options=optim_options)
        opt_solution.update({"sample_name": [initial_condition_id] * len(opt_solution["data_set_id"])})
        # update dictionary
        for key, value in it.chain(opt_solution.items(), empty_dict.items()):
            for i_value in value:
                all_initial_condition_sol[key].append(i_value)
        # process opt solution for ach initial condition for all experimental data sets
        # all_parameter_info = process_opt_solution(opt_solution, number_of_parameters=number_of_parameters,
        #                                           flux_id=flux_id, flux_choice=flux_choice)
        # all_x0_all_parameter_opt_info.append(all_parameter_info)
        print("Analysis with Initial Condition {} of {} Complete".format(j_id+1, number_initial_conditions))

    # process and write data to file
    all_sol_df, index_labels = write_ident_info_file(all_initial_condition_sol, exp_df, file_name)

    return all_sol_df, index_labels


def generate_random_initial_conditions(given_initial_condition, number_random_conditions, negative=0):
    """generate random initial conditions based around given initial conidtion"""
    number_variables = len(given_initial_condition)
    # set random number generator seed
    rnd_num = RandomState(12345678)
    if negative:
        pos_rnd_value_changes = rnd_num.uniform(0, 1, (number_variables, number_random_conditions/2))
        neg_rnd_value_changes = - rnd_num.uniform(0, 1, (number_variables, number_random_conditions / 2))
        rnd_value_changes = np.hstack((pos_rnd_value_changes, neg_rnd_value_changes))
    else:
        rnd_value_changes = rnd_num.uniform(0, 1, (number_variables, number_random_conditions))

    given_initial_condition = np.reshape(given_initial_condition, (6, 1))
    new_initial_conditions = np.repeat(given_initial_condition, 10, axis=1) * (1 + rnd_value_changes)

    new_initial_conditions = list(np.transpose(new_initial_conditions))
    return new_initial_conditions
