from mpi4py import MPI
from mpi_master_slave import Master, Slave
from mpi_master_slave import WorkQueue
import numpy as np


class ParallelFlux(object):
    """Class running multiple ode simulations by assigning jobs to different slaves
    until all the work is done"""

    def __init__(self, slaves):
        # when creating the Master we tell it what slaves it can handle
        self.master = Master(slaves)
        # WorkQueue is a convenient class that run slaves on a tasks queue
        self.work_queue = WorkQueue(self.master)

    def terminate_slaves(self):
        """
        Call this to make all slaves exit their run loop
        """
        self.master.terminate_slaves()

    def __add_next_task(self, flux_fun, concentration, parameters, value_id=()):
        """
        create tasks and add it to the work queue
        Every task has specific arguments
        """
        # set ode rhs function, initial condition and parameters
        data = {'flux_fun': flux_fun, 'y': concentration, 'id': value_id, 'ode_sys_opts': parameters}

        # add data to work queue
        self.work_queue.add_work(data)

    def run_i_value(self, flux_fun, concentration, ode_sys_opts):
        """
        This is the core where I keep starting slaves
        as long as there is work to do
        """
        # let's prepare our work queue. This can be built at initialization time
        # but it can also be added later as more work become available
        #
        for idx, (j_conc_value, j_parameter) in enumerate(zip(concentration, ode_sys_opts)):
            self.__add_next_task(flux_fun=flux_fun, concentration=j_conc_value, parameters=j_parameter, value_id=idx)

        # Keep starting slaves as long as there is work to do
        all_boolean = []
        all_flux = []
        all_id = []
        while not self.work_queue.done():
            # give more work to do to each idle slave (if any)
            self.work_queue.do_work()
            # reclaim returned data from completed slaves
            for slave_return_data in self.work_queue.get_completed_work():
                y0_id, done, j_slave_flux = slave_return_data
                # import pdb; pdb.set_trace()
                all_boolean.append(done)
                all_flux.append(j_slave_flux)
                all_id.append(y0_id)

                if done:
                    print('Master: slave finished its task returning: %s)' % str(y0_id))
            # sleep some time
            # time.sleep(0.3)
        results = {'flux': all_flux, 'y0_id': all_id, 'boolean': all_boolean}
        return results


class MySlave(Slave):
    """
    A slave process extends Slave class, overrides the 'do_work' method
    and calls 'Slave.run'. The Master will do the rest
    In this example we have different tasks but instead of creating a Slave for
    each type of taks we create only one class that can handle any type of work.
    This avoids having idle processes if, at certain times of the execution, there
    is only a particular type of work to do but the Master doesn't have the right
    slave for that task.
    """

    def __init__(self):
        super(MySlave, self).__init__()

    def do_work(self, data):
        """do work method overrides Slave.do_work() and defines work
        to be done by every slave"""
        # boolean variable set to determine whether work on slave was completed on not
        done = False

        # define explicit assimulo problem
        flux_fun = data['flux_fun']
        value_id = data['id']
        ode_sys_opts = data['ode_sys_opts']
        concentration = data['concentration']
        # flux_value = flux_fun(concentration, ode_sys_opts)
        flux_value = np.array(list(map(lambda x: flux_fun(x, ode_sys_opts), concentration)))

        rank = MPI.COMM_WORLD.Get_rank()
        name = MPI.Get_processor_name()

        print('  Slave %s rank %d executing task %s' % (name, rank, value_id))

        if len(flux_value):
            done = True

        if done:
            return value_id, done, flux_value
        else:
            return value_id, done, []


def setup_parallel_flux(flux_fun, concentrations, parameters):

    name = MPI.Get_processor_name()
    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()

    all_flux = {}
    print('I am  %s rank %d (total %d)' % (name, rank, size))
    if rank == 0:  # Master
        ode_job = ParallelFlux(slaves=range(1, size))
        all_flux = ode_job.run_i_value(flux_fun=flux_fun, concentration=concentrations, ode_sys_opts=parameters)
        ode_job.terminate_slaves()
    else:  # Any slave
        MySlave().run()
    import pdb; pdb.set_trace()
    return all_flux
