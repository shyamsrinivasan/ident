import numpy as np
from assimulo.solvers import CVode
from assimulo.problem import Explicit_Problem


def simulate_ode(fun, y_initial, tf, opts):
    "function to run CVode solver on given problem"
    # get options
    iter, discretization_method, atol, rtol, time_points, ode_system_options = opts

    ode_function = lambda t, x : fun(t,x,ode_system_options)

    # define explicit assimulo problem
    prob = Explicit_Problem(ode_function, y0=y_initial)

    # create solver instance
    solver = CVode(prob)

    # set solver options
    solver.iter, solver.discr, solver.atol, solver.rtol = iter, discretization_method, atol, rtol

    # simulate system
    time_course, y_result = solver.simulate(tf, time_points)

    return time_course, y_result, prob, solver
