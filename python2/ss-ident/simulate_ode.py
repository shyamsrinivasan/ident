import numpy as np
import matplotlib.pyplot as plt
from assimulo.solvers import CVode
from assimulo.problem import Explicit_Problem


def simulate_ode(fun, y_initial, tf, opts):
    "function to run CVode solver on given problem"
    # get options
    ode_opts, ode_system_options = opts
    iter, discretization_method, atol, rtol, time_points = ode_opts

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

def run_ode_sims(fun, y_initial, opts, t_final=500, args_1=False):
    """run kotte model ode using cvode from assimulo"""

    # def_par_val = np.array([.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, .1])
    time_points, y_dynamic, prob, solver = simulate_ode(fun, y_initial, t_final, opts)
    if args_1:
        plt.plot(time_points, y_dynamic, color="r")
        plt.xlabel('Time')
        plt.ylabel('Dependent Variables')
        plt.show()

    return time_points, y_dynamic, prob, solver
