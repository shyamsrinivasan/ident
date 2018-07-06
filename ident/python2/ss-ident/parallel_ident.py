from mpi4py import MPI
from mpi_master_slave import Master, Slave
from mpi_master_slave import WorkQueue
from identifiability_analysis import run_flux_ident
import pandas as pd
# import numpy as np
import time


class ParallelIdent(object):
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

    def __add_next_task(self, exp_data, ident_fun, flux_id, flux_choice, sample_id, data_set_id):
        """
        create tasks and add it to the work queue
        Every task has specific arguments
        """
        # set ode rhs function, initial condition and parameters
        data = {'ident_fun': ident_fun, 'flux_id': flux_id, 'flux_choice': flux_choice,
                'sample_id': sample_id, 'data_set_id': data_set_id, 'exp_data': exp_data}
        # add data to work queue
        self.work_queue.add_work(data)

    def run_i_value(self, exp_df, ident_fun, flux_id, flux_choice):
        """
        This is the core where I keep starting slaves
        as long as there is work to do
        """

        # let's prepare our work queue. This can be built at initialization time
        # but it can also be added later as more work become available
        #
        idx = pd.IndexSlice
        all_df_indices = exp_df.index.unique().tolist()
        # create tuple of indices
        # import pdb; pdb.set_trace()
        for j_index, sample_data_set_id in enumerate(all_df_indices):
            j_exp_data_set = exp_df.loc[idx[sample_data_set_id],
                                        ['acetate', 'pep', 'fdp', 'E', 'v1', 'v2', 'v3', 'v5']].values.tolist()
            flat_data_list = [i_element for i_data in j_exp_data_set for i_element in i_data]
            self.__add_next_task(exp_data=flat_data_list, ident_fun=ident_fun, flux_id=flux_id, flux_choice=flux_choice,
                                 sample_id=sample_data_set_id[0], data_set_id=sample_data_set_id[1])

        # Keeep starting slaves as long as there is work to do
        ident_results = []
        while not self.work_queue.done():
            # give more work to do to each idle slave (if any)
            self.work_queue.do_work()
            # reclaim returned data from completed slaves
            for slave_return_data in self.work_queue.get_completed_work():
                complete, ident_info, slave_flux_id, slave_flux_choice, sample_id, data_set_id = slave_return_data
                i_slave_result = {'done': complete, 'flux_id': flux_id, 'flux_choice': flux_choice, 'sample_id': sample_id,
                                  'data_set_id': data_set_id, 'ident_info': ident_info}
                ident_results.append(i_slave_result)

                if complete:
                    print('Master: slave finished its task returning: %s)' % str(data_set_id))
            # sleep some time
            # time.sleep(0.01)
        return ident_results


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

        # define explicit assimulo problem
        ident_fun = data['ident_fun']
        exp_data = data['exp_data']
        flux_id = data['flux_id']
        flux_choice = data['flux_choice']
        sample_id = data['sample_id']
        data_set_id = data['data_set_id']
        try:
            ident_info, _, _ = run_flux_ident(ident_fun, exp_data, flux_id, flux_choice)
            complete = True
        except TypeError:
            ident_info = []
            complete = False

        rank = MPI.COMM_WORLD.Get_rank()
        name = MPI.Get_processor_name()

        print('  Slave %s rank %d executing task (sample: %s, data set: %s)' % (name, rank, sample_id, data_set_id))

        if complete:
            return complete, ident_info, flux_id, flux_choice, sample_id, data_set_id
        else:
            return complete, [], [], [], [], []


def setup_parallel_ident(ident_fun, flux_id, flux_choice, exp_data):
    name = MPI.Get_processor_name()
    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()

    sim_result = {}
    print('I am  %s rank %d (total %d)' % (name, rank, size))
    if rank == 0:  # Master
        ident_job = ParallelIdent(slaves=range(1, size))
        # import pdb;pdb.set_trace()
        sim_result = ident_job.run_i_value(exp_data, ident_fun=ident_fun, flux_id=flux_id, flux_choice=flux_choice)

        ident_job.terminate_slaves()
    else:  # Any slave
        MySlave().run()

    return sim_result
