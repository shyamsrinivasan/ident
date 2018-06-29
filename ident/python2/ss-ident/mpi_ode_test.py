import numpy as np
from parallel_ode import ParallelOde


def py_rhs_fun(t, y, p):
    alpha = p[0]
    beta = p[1]
    delta = p[2]
    gamma = p[3]

    flux = np.array([alpha * y[0], beta * y[0] * y[1], delta * y[0] * y[1], gamma * y[1]])
    rhs = np.array([flux[0] - flux[1], flux[2] - flux[3]])

    return rhs


if __name__ == "__main__":
    name = MPI.Get_processor_name()
    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()

    print('I am  %s rank %d (total %d)' % (name, rank, size))

    if rank == 0:  # Master
        # import pdb; pdb.set_trace()
        # set parameter value
        parameters = {'ode_opts': {'iter': 'Newton', 'discr': 'Adams', 'atol': 1e-10, 'rtol': 1e-10,
                                   'time_points': 200, 'display_progress': True, 'verbosity': 30},
                      'ode_sys_opts': np.array([.5, .02, .4, .004])}
        y0 = [np.array([1, .0001]), np.array([2, .0001]), np.array([3, .0001]), np.array([4, .0001]),
              np.array([5, .0001]), np.array([.0001, 1]), np.array([.0001, 2]), np.array([.0001, 3]),
              np.array([.0001, 4]), np.array([.0001, 5]), np.array([10, 1]), np.array([1, 10])]
        # import pdb; pdb.set_trace()
        ode_job = ParallelOde(slaves=range(1, size))
        all_boolean, all_tout, all_yout, all_y0_id = ode_job.run(ode_rhs_fun=py_rhs_fun, t_final=200,
                                                                 parameters=parameters, initial_values=y0)
        import pdb; pdb.set_trace()
        app.terminate_slaves()

    else:  # Any slave

        MySlave().run()

    print('Task completed (rank %d)' % (rank))