import numpy as np
from assimulo.solvers import CVode
from assimulo.problem import Explicit_Problem


def simulate_ode(fun, y_initial, tf, opts):
    "function to run CVode solver on given problem"
    # define explicit assimulo problem
    prob = Explicit_Problem(fun, y0=y_initial)

    # create solver instance
    solver = CVode(prob)

    # set solver options
    solver.iter, solver.discr, solver.atol, solver.rtol, time_points = opts

    # simulate system
    time_course, y_result = solver.simulate(tf, time_points)

    return time_course, y_result, prob, solver