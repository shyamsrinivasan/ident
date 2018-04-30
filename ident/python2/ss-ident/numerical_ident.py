import casadi as casadi


def v3_ck_numerical_problem(chosen_data):
    """fun to calculate objectives and constraints for numerical identification of v3
    using convenience kinetics"""

    # symbolic variables for use with casadi
    v3max = casadi.SX.sym("v3max")
    k3fdp = casadi.SX.sym("k3fdp")
    k3pep = casadi.SX.sym("k3pep")
    error = casadi.SX.sym("error", 3)

    all_constraints = []
    for i_data in range(0, 3):
        _, pep, fdp, _, _, _, v3, _, _, _ = list(chosen_data[i_data])

        # flux equation for each experiment i_data in chosen_data
        fdp_sat = fdp / k3fdp
        pep_sat = pep / k3pep
        nr_3 = v3max * fdp_sat
        dr_3 = 1 + fdp_sat
        regulation_activate = 1 / (1 + 1 / pep_sat)
        flux_3 = regulation_activate * nr_3 / dr_3
        # formulate constraint for each experimental data set
        all_constraints.append(flux_3 - v3 - error[i_data])

    # create casadi array of constraints
    constraints = casadi.vertcat(*all_constraints)

    # create complete casadi objective fun
    objective = error[0]**2
    for i_error_id in range(1, 3):
        objective += error[i_error_id]**2

    # get complete set of all optimization variables vertcat([parameters, error])
    x = casadi.vertcat(*[v3max, k3fdp, k3pep, error])

    nlp = {'x': x, 'f': objective, 'g': constraints}

    return nlp


def v3_numerical_problem(chosen_data):
    """fun to calculate objectives and constraints for numerical identification of v3
        using MWC kinetics"""

    # symbolic variables for use with casadi
    v3max = casadi.SX.sym("v3max")
    k3fdp = casadi.SX.sym("k3fdp")
    k3pep = casadi.SX.sym("k3pep")
    l3fdp = casadi.SX.sym("l3fdp")
    error = casadi.SX.sym("error", 4)

    # fdp_sat = 1 + y[1] / K3fdp
    # pep_sat = 1 + y[0] / K3pep
    # flux_3 = V3max * (fdp_sat - 1) * (fdp_sat ** 3) / (fdp_sat ** 4 + L3fdp * (pep_sat ** (-4)))

    all_constraints = []
    for i_data in range(0, 4):
        _, pep, fdp, _, _, _, v3, _, _, _ = list(chosen_data[i_data])

        # flux equation for each experiment i_data in chosen_data
        fdp_sat = 1 + fdp / k3fdp
        pep_sat = 1 + pep / k3pep
        nr_3 = v3max * (fdp_sat - 1) * (fdp_sat**3)
        dr_3 = (fdp_sat ** 4 + l3fdp * 1e6 * pep_sat ** (-4))
        flux_3 = nr_3 / dr_3
        # formulate constraint for each experimental data set
        all_constraints.append(flux_3 - v3 - error[i_data])

    # create casadi array of constraints
    constraints = casadi.vertcat(*all_constraints)

    # create complete casadi objective fun
    objective = error[0] ** 2
    for i_error_id in range(1, 4):
        objective += error[i_error_id] ** 2

    # get complete set of all optimization variables vertcat([parameters, error])
    x = casadi.vertcat(*[v3max, k3fdp, k3pep, l3fdp, error])

    nlp = {'x': x, 'f': objective, 'g': constraints}

    return nlp


def solve_numerical_nlp(chosen_fun, chosen_data):
    """solve nlp for determining parameter values numerically"""
    # formulate nlp
    arg = {}
    if chosen_fun == 0:
        nlp = v3_ck_numerical_problem(chosen_data)
        # Initial condition
        arg["x0"] = [.1, .1, .1, 0, 0, 0]
        # Bounds on x
        ubx = [2, 1, 1, .1, .1, .1]
        # Bounds on the constraints
        arg["lbg"] = 3 * [0]
        arg["ubg"] = 3 * [0]
    elif chosen_fun == 1:
        nlp = v3_numerical_problem(chosen_data)
        # Initial condition
        arg["x0"] = [.1, .1, .1, 1, 0, 0, 0, 0]
        # Bounds on x
        ubx = [2, 1, 1, 5, .1, .1, .1, .1]
        # Bounds on the constraints
        arg["lbg"] = 4 * [0]
        arg["ubg"] = 4 * [0]

    # NLP solver options
    opts = {"ipopt.tol": 1e-12}
    # opts = {"qpsol": "qpoases"}

    # Allocate an NLP solver and buffer
    solver = casadi.nlpsol("solver", "ipopt", nlp, opts)
    # solver = casadi.nlpsol("solver", "sqpmethod", nlp, opts)

    # Bounds on x
    lbx = nlp["x"].nnz() * [0]

    arg["lbx"] = lbx
    arg["ubx"] = ubx

    # Solve the problem
    res = solver(**arg)

    # Print the optimal cost
    print("optimal cost: ", float(res["f"]))

    # Print the optimal solution
    xopt = res["x"]
    print("optimal solution: ", xopt)
    print("dual solution (x) = ", res["lam_x"])
    print("dual solution (g) = ", res["lam_g"])

    return None
