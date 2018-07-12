from mpi4py import MPI
from mpi_master_slave import Master, Slave
from mpi_master_slave import WorkQueue
from simulate_ode import simulate_ode
from names_strings import true_parameter_values
from run_validate import ValidateSim
from run_ident import ModelIdent
import numpy as np
import kotte_model
import os.path
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
        """,
        Call this to make all slaves exit their run loop
        """
        self.master.terminate_slaves()

    def __add_next_task(self, task, **kwargs):
        """
        create tasks and add it to the work queue. Every task has specific arguments
        """

        if task == 'initial_sim':
            data = {'sim_obj': kwargs['sim_obj'], 'ode_sys_opts': kwargs['ode_sys_opts'],
                    'id': kwargs['estimate_id'], 'y0': kwargs['y0'], 'task': task}
            # 'ode_fun': kwargs['rhs_fun'],
            # 'ode_opts': kwargs['ode_opts'],
            # 't_final': kwargs['t_final'], 'flux_fun': kwargs['flux_fun'], }

        self.work_queue.add_work(data)

    def run_all(self, task, **kwargs):
        """parallel application core"""

        if task == 'initial_sim':
            import pdb; pdb.set_trace()
            estimates = kwargs['parameters']
            estimate_id = kwargs['estimate_info']
            sim_obj = kwargs['sim_obj']
            for j_index, (j_estimate, j_estimate_id) in enumerate(zip(estimates, estimate_id)):
                self.__add_next_task(task=task, **{'sim_obj': sim_obj, 'ode_sys_opts': j_estimate,
                                                   'estimate_id': j_estimate_id, 'y0': sim_obj.wt_y0})
                                                   # 'rhs_fun': sim_obj.rhs_fun,
                                                   # 'flux_fun': sim_obj.flux_fun, 't_final': sim_obj.t_final,
                                                   #  'ode_opts': sim_obj.ode_opts})

        # Keeep starting slaves as long as there is work to do
        results = []
        while not self.work_queue.done():
            # give more work to do to each idle slave (if any)
            self.work_queue.do_work()

            # reclaim returned data from completed slaves
            for slave_return_data in self.work_queue.get_completed_work():
                task, output = slave_return_data
                if task == 'initial_sim':
                    j_slave_t, j_slave_y, j_slave_f, estimate_id, sample_id, data_set_id = output
                    j_slave_dyn = {'t': j_slave_t, 'y': j_slave_y, 'flux': j_slave_f, 'estimate_id': estimate_id,
                                   'sample_id': sample_id, 'data_set_id': data_set_id}

                    # get ss values
                    j_slave_ss = {'y': j_slave_y[-1, :], 'flux': j_slave_f[-1, :], 'estimate_id': estimate_id,
                                  'sample_id': sample_id, 'data_set_id': data_set_id}
                    i_slave_result = {'dynamic': j_slave_dyn, 'ss': j_slave_ss}
                    results.append(i_slave_result)

        return results


class ValidateSlave(Slave):
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
        super(ValidateSlave, self).__init__()

    def do_work(self, data):
        """do work method overrides Slave.do_work() and defines work
        to be done by every slave"""
        rank = MPI.COMM_WORLD.Get_rank()
        name = MPI.Get_processor_name()

        print('  Slave %s rank %d executing task %s' % (name, rank, data['task']))

        if data['task'] == 'initial_sim':
            # define explicit assimulo problem
            sim_obj = data['sim_obj']
            rhs_fun = sim_obj.rhs_fun  # data['rhs_fun']
            y_initial = data['y0']
            estimate_id = data['estimate_id']
            ode_opts = sim_obj.ode_opts  # data['ode_opts']
            ode_sys_opts = data['ode_sys_opts']
            t_final = sim_obj.t_final  # data['t_final']
            all_options = [ode_opts, ode_sys_opts]

            print('  Slave %s rank %d executing initial_sim for estimate: %s sample: %s, data set: %s' %
                  (name, rank, estimate_id[0], estimate_id[1], estimate_id[2]))
            slave_tout, slave_yout, _, _ = simulate_ode(rhs_fun, y_initial, tf=t_final, opts=all_options)
            print(' ode simulation complete ')

            # calculate flux
            flux_fun = sim_obj.flux_fun  # data['flux_fun']
            slave_flux = np.array(list(map(lambda x: flux_fun(x, ode_sys_opts), slave_yout)))

            return data['task'], (slave_tout, slave_yout, slave_flux, estimate_id[0], estimate_id[1], estimate_id[2])


# class MySlave(Slave):
#     """
#     A slave process extends Slave class, overrides the 'do_work' method
#     and calls 'Slave.run'. The Master will do the rest
#     In this example we have different tasks but instead of creating a Slave for
#     each type of taks we create only one class that can handle any type of work.
#     This avoids having idle processes if, at certain times of the execution, there
#     is only a particular type of work to do but the Master doesn't have the right
#     slave for that task.
#     """
#
#     def __init__(self):
#         super(MySlave, self).__init__()
#
#     @staticmethod
#     def change_parameter_values(changed_parameter, default_parameter):
#         """change default parameter by value metnioned in one of changed parameter"""
#         # if not default_parameter:
#         #     default_parameter = self.i_parameter
#
#         new_parameter_list = []
#         for i_parameter_change in changed_parameter:
#             new_parameter = deepcopy(default_parameter)
#             parameter_name = list(i_parameter_change.keys())[0]
#             parameter_change = np.array(list(i_parameter_change.values())[0])
#             if parameter_name == 'wt':
#                 new_parameter['ac'] = new_parameter['ac'] * (1 + parameter_change)
#             else:
#                 new_parameter[parameter_name] = new_parameter[parameter_name] * (1 + parameter_change)
#             new_parameter_list.append(new_parameter)
#         return new_parameter_list
#
#     def do_work(self, data):
#         """do work method overrides Slave.do_work() and defines work
#         to be done by every slave"""
#
#         dyn_results = []
#         ss_results = []
#
#         # do initial simulation run
#         rhs_fun = data['ode_fun']
#         y_initial = data['y0']
#         estimate_id = data['id']
#         ode_opts = data['ode_opts']
#         ode_sys_opts = data['ode_sys_opts']
#         t_final = data['t_final']
#         all_options = [ode_opts, ode_sys_opts]
#         tout_i, yout_i, _, _ = simulate_ode(rhs_fun, y_initial, tf=t_final, opts=all_options)
#
#         # calculate flux
#         flux_fun = data['flux_fun']
#         fout_i = np.array(list(map(lambda x: flux_fun(x, ode_sys_opts), yout_i)))
#
#         # get ss values
#         ss_results.append({'y': yout_i[-1, :], 'flux': fout_i[-1, :], 'id': 'initial_ss'})
#         dyn_results.append({'y': yout_i, 'flux': fout_i, 't': tout_i, 'id': 'initial_ss'})
#
#         # create list of perturbed parameters
#         new_parameter_list = self.change_parameter_values(changed_parameter=data['perturbations'],
#                                                           default_parameter=ode_sys_opts)
#
#         # do perturbation simulation run
#         for exp_id, i_parameter in enumerate(new_parameter_list):
#             print(' Perturbation %d \n' % exp_id)
#             all_options = [ode_opts, i_parameter]
#             tout_p, yout_p, _, _ = simulate_ode(rhs_fun, yout_i[-1, :], tf=t_final, opts=all_options)
#             fout_p = np.array(list(map(lambda x: flux_fun(x, i_parameter), yout_p)))
#
#             # collect results
#             ss_results.append({'y': yout_p[-1, :], 'flux': fout_p[-1, :], 'id': 'experiment_{}'.format(exp_id)})
#             dyn_results.append({'y': yout_p, 'flux': fout_p, 't': tout_p, 'id': 'experiment_{}'.format(exp_id)})
#
#         rank = MPI.COMM_WORLD.Get_rank()
#         name = MPI.Get_processor_name()
#
#         print('  Slave %s rank %d executing task %s' % (name, rank, estimate_id))
#
#         return estimate_id, ss_results, dyn_results


def v1_validate():
    """validate estimated parameter values"""

    name = MPI.Get_processor_name()
    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()

    if rank == 0:
        # create ident object first
        v1_ident = ModelIdent(ident_fun=kotte_model.flux_1_kcat_ident,
                              arranged_data_file_name=os.path.join(os.getcwd(), 'exp/exp_v1_2_experiments'),
                              ident_data_file_name=os.path.join(os.getcwd(), 'ident/ident_v1_kcat'),
                              **{'original_exp_file': os.path.join(os.getcwd(), 'exp/experiments'),
                                 'flux_id': 1, 'flux_choice': 2,
                                 'values_figure': os.path.join(os.getcwd(), 'results/v1_kcat_parameter_values.eps'),
                                 'ident_figure': os.path.join(os.getcwd(), 'results/v1_kcat_ident.eps'),
                                 'exp_figure': os.path.join(os.getcwd(), 'results/v1_kcat_exp.eps'),
                                 'figure_format': 'eps'})

        # retrieve identifiability data and process it for validation
        v1_ident.validation_info()

        user_ode_opts = {'iter': 'Newton', 'discr': 'Adams', 'atol': 1e-10, 'rtol': 1e-10,
                         'time_points': 200, 'display_progress': True, 'verbosity': 30}
        # initial ss to begin all simulations from
        y0 = np.array([5, 1, 1])
        # get and set true parameter values, if available separately
        default_parameters = true_parameter_values()

        v1_validate = ValidateSim(kotte_model.kotte_ck_ode, kotte_model.kotte_ck_flux, **{'kinetics': 2,
                                                                                          'ode_opts': user_ode_opts,
                                                                                          't_final': 200,
                                                                                          'wt_y0': y0,
                                                                                          'i_parameter':
                                                                                              default_parameters,
                                                                                          'sample_size': 1,
                                                                                          'noise_std': 0.05})

        parameter_estimates, estimate_info = v1_validate.create_parameter_list(v1_ident.select_values)

        job = ParallelValidate(slaves=range(1, size))

        ident_result = job.run_all(task='initial_sim', **{'parameters': parameter_estimates,
                                                          'estimate_info': estimate_info, 'sim_obj': v1_validate})
        job.terminate_slaves()

    else:
        print('I am %s Slave with rank %s of %s' % (name, str(rank), str(size)))
        ValidateSlave().run()



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
