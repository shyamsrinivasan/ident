from mpi4py import MPI
from mpi_master_slave import Master, Slave
from mpi_master_slave import WorkQueue
from simulate_ode import simulate_ode
import numpy as np
from copy import deepcopy


class ParallelValidate(object):
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

    def __add_next_task(self, j_validate_object, perturbations):
        """
        create tasks and add it to the work queue
        Every task has specific arguments
        """
        # set ode rhs function, initial condition and parameters
        data = {'ode_fun': j_validate_object.ode_rhs_fun, 'y0': j_validate_object.initial_value,
                'id': j_validate_object.estimate_id, 'ode_sys_opts': j_validate_object.i_parameter,
                'ode_opts': j_validate_object.ode_opts, 't_final': j_validate_object.t_final,
                'flux_fun': j_validate_object.flux_fun, 'perturbations': perturbations}

        # add data to work queue
        self.work_queue.add_work(data)

    def run_i_value(self, validate_sim_objs, perturbations):
        """
        This is the core where I keep starting slaves
        as long as there is work to do
        """

        # let's prepare our work queue. This can be built at initialization time
        # but it can also be added later as more work become available
        #
        for j_object in validate_sim_objs:
            self.__add_next_task(j_validate_object=j_object, perturbations=perturbations)

        # Keeep starting slaves as long as there is work to do
        all_id = []
        all_ss = []
        all_dyn = []
        while not self.work_queue.done():
            # give more work to do to each idle slave (if any)
            self.work_queue.do_work()
            # reclaim returned data from completed slaves
            for slave_return_data in self.work_queue.get_completed_work():
                estimate_id, ss_results, dyn_results = slave_return_data
                # import pdb; pdb.set_trace()
                all_id.append(estimate_id)
                all_ss.append(ss_results)
                all_dyn.append(dyn_results)

                print('Master: slave finished its task returning: %s)' % str(estimate_id))
            # sleep some time
            # time.sleep(0.3)
        results = {'estimate_id': all_id, 'ss_info': all_ss, 'dyn_info': all_dyn}
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

    @staticmethod
    def change_parameter_values(changed_parameter, default_parameter):
        """change default parameter by value metnioned in one of changed parameter"""
        # if not default_parameter:
        #     default_parameter = self.i_parameter

        new_parameter_list = []
        for i_parameter_change in changed_parameter:
            new_parameter = deepcopy(default_parameter)
            parameter_name = list(i_parameter_change.keys())[0]
            parameter_change = np.array(list(i_parameter_change.values())[0])
            if parameter_name == 'wt':
                new_parameter['ac'] = new_parameter['ac'] * (1 + parameter_change)
            else:
                new_parameter[parameter_name] = new_parameter[parameter_name] * (1 + parameter_change)
            new_parameter_list.append(new_parameter)
        return new_parameter_list

    def do_work(self, data):
        """do work method overrides Slave.do_work() and defines work
        to be done by every slave"""

        dyn_results = []
        ss_results = []

        # do initial simulation run
        rhs_fun = data['ode_fun']
        y_initial = data['y0']
        estimate_id = data['id']
        ode_opts = data['ode_opts']
        ode_sys_opts = data['ode_sys_opts']
        t_final = data['t_final']
        all_options = [ode_opts, ode_sys_opts]
        tout_i, yout_i, _, _ = simulate_ode(rhs_fun, y_initial, tf=t_final, opts=all_options)

        # calculate flux
        flux_fun = data['flux_fun']
        fout_i = np.array(list(map(lambda x: flux_fun(x, ode_sys_opts), yout_i)))

        # get ss values
        ss_results.append({'y': yout_i[-1, :], 'flux': fout_i[-1, :], 'id': 'initial_ss'})
        dyn_results.append({'y': yout_i, 'flux': fout_i, 't': tout_i, 'id': 'initial_ss'})

        # create list of perturbed parameters
        new_parameter_list = self.change_parameter_values(changed_parameter=data['perturbations'],
                                                          default_parameter=ode_sys_opts)

        # do perturbation simulation run
        for exp_id, i_parameter in enumerate(new_parameter_list):
            print(' Perturbation %d \n' % exp_id)
            all_options = [ode_opts, i_parameter]
            tout_p, yout_p, _, _ = simulate_ode(rhs_fun, yout_i[-1, :], tf=t_final, opts=all_options)
            fout_p = np.array(list(map(lambda x: flux_fun(x, i_parameter), yout_p)))

            # collect results
            ss_results.append({'y': yout_p[-1, :], 'flux': fout_p[-1, :], 'id': 'experiment_{}'.format(exp_id)})
            dyn_results.append({'y': yout_p, 'flux': fout_p, 't': tout_p, 'id': 'experiment_{}'.format(exp_id)})

        rank = MPI.COMM_WORLD.Get_rank()
        name = MPI.Get_processor_name()

        print('  Slave %s rank %d executing task %s' % (name, rank, estimate_id))

        return estimate_id, ss_results, dyn_results


def setup_parallel_validate(validate_sim_objects, perturbations):

    name = MPI.Get_processor_name()
    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()

    sim_result = {}
    print('I am  %s rank %d (total %d)' % (name, rank, size))
    if rank == 0:  # Master
        validate_job = ParallelValidate(slaves=range(1, size))
        # import pdb;pdb.set_trace()
        sim_result = validate_job.run_i_value(validate_sim_objects, perturbations)

        validate_job.terminate_slaves()
    else:  # Any slave
        MySlave().run()

    # rearrange results based on order of parameters/initial values
    # import pdb;pdb.set_trace()
    return sim_result
