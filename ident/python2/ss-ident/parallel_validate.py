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

        elif task == 'perturbation_sim':
            data = {'sim_obj': kwargs['sim_obj'], 'ode_sys_opts': kwargs['ode_sys_opts'],
                    'id': kwargs['estimate_id'], 'y0': kwargs['y0'], 'perturbation_id': kwargs['perturbation_id'],
                    'task': task}

        self.work_queue.add_work(data)

    def run_all(self, task, **kwargs):
        """parallel application core"""

        if task == 'initial_sim':
            estimates = kwargs['parameters']
            estimate_id = kwargs['estimate_info']
            sim_obj = kwargs['sim_obj']
            for j_index, (j_estimate, j_estimate_id) in enumerate(zip(estimates, estimate_id)):
                self.__add_next_task(task=task, **{'sim_obj': sim_obj, 'ode_sys_opts': j_estimate,
                                                   'estimate_id': j_estimate_id, 'y0': sim_obj.wt_y0})

        # Keeep starting slaves as long as there is work to do
        results = []
        while not self.work_queue.done():
            # give more work to do to each idle slave (if any)
            self.work_queue.do_work()

            # reclaim returned data from completed slaves
            for slave_return_data in self.work_queue.get_completed_work():
                task, output = slave_return_data

                if task == 'initial_sim':
                    j_slave_t, j_slave_y, j_slave_f, estimate_id, sample_id, data_set_id, sim_obj, ode_sys_opts = output
                    j_slave_dyn = {'t': j_slave_t, 'y': j_slave_y, 'flux': j_slave_f, 'estimate_id': estimate_id,
                                   'sample_id': sample_id, 'data_set_id': data_set_id}

                    # info on bistability
                    if j_slave_y[-1, 0] > j_slave_y[-1, 1]:
                        stability_id = 1
                    elif j_slave_y[-1, 0] < j_slave_y[-1, 1]:
                        stability_id = 2
                    else:
                        stability_id = 0

                    # get ss values
                    j_slave_ss = {'y': j_slave_y[-1, :], 'flux': j_slave_f[-1, :], 'estimate_id': estimate_id,
                                  'sample_id': sample_id, 'data_set_id': data_set_id, 'ssid': stability_id}
                    i_slave_result = {'dynamic': j_slave_dyn, 'ss': j_slave_ss, 'initial': True, 'perturbation': False}

                    # append initial sim results
                    results.append(i_slave_result)

                    # create list of perturbated parameter for given parameter estimate in ode_sys_opts
                    perturbed_parameter_list = sim_obj.change_parameter_values(sim_obj.test_perturbations, ode_sys_opts)

                    # add perturbation jobs to back of jobs queue
                    experiment_id = ['experiment_{}'.format(parameter_id) for parameter_id, _ in
                                     enumerate(perturbed_parameter_list)]
                    for j_experiment_id, j_perturbation in zip(experiment_id, perturbed_parameter_list):
                        self.__add_next_task(task='perturbation_sim', **{'sim_obj': sim_obj,
                                                                         'ode_sys_opts': j_perturbation, 'estimate_id':
                                                                             (estimate_id, sample_id, data_set_id),
                                                                         'perturbation_id': j_experiment_id,
                                                                         'y0': j_slave_y[-1, :]})

                elif task == 'perturbation_sim':
                    j_slave_t, j_slave_y, j_slave_f, estimate_id, sample_id, data_set_id, perturbation_id = output
                    j_slave_dyn = {'t': j_slave_t, 'y': j_slave_y, 'flux': j_slave_f, 'estimate_id': estimate_id,
                                   'sample_id': sample_id, 'data_set_id': data_set_id, 'perturbation_id':
                                       perturbation_id}

                    # info on bistability
                    if j_slave_y[-1, 0] > j_slave_y[-1, 1]:
                        stability_id = 1
                    elif j_slave_y[-1, 0] < j_slave_y[-1, 1]:
                        stability_id = 2
                    else:
                        stability_id = 0

                    # get ss values
                    j_slave_ss = {'y': j_slave_y[-1, :], 'flux': j_slave_f[-1, :], 'estimate_id': estimate_id,
                                  'sample_id': sample_id, 'data_set_id': data_set_id, 'perturbation_id':
                                      perturbation_id, 'ssid': stability_id}
                    i_slave_result = {'dynamic': j_slave_dyn, 'ss': j_slave_ss, 'initial': False, 'perturbation': True}

                    # append perturbation results
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
            estimate_id = data['id']
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

            result = (slave_tout, slave_yout, slave_flux, estimate_id[0], estimate_id[1], estimate_id[2], sim_obj,
                      ode_sys_opts)

        elif data['task'] == 'perturbation_sim':

            sim_obj = data['sim_obj']
            rhs_fun = sim_obj.rhs_fun  # data['rhs_fun']
            y_initial = data['y0']
            estimate_id = data['id']
            perturbation_id = data['perturbation_id']
            ode_opts = sim_obj.ode_opts  # data['ode_opts']
            ode_sys_opts = data['ode_sys_opts']
            t_final = sim_obj.t_final  # data['t_final']
            all_options = [ode_opts, ode_sys_opts]

            print('  Slave %s rank %d executing initial_sim for estimate: %s sample: %s, data set: %s '
                  'perturbation: %s' %
                  (name, rank, estimate_id[0], estimate_id[1], estimate_id[2], perturbation_id))
            slave_tout, slave_yout, _, _ = simulate_ode(rhs_fun, y_initial, tf=t_final, opts=all_options)
            print(' ode perturbation simulation complete ')

            # calculate flux
            flux_fun = sim_obj.flux_fun  # data['flux_fun']
            slave_flux = np.array(list(map(lambda x: flux_fun(x, ode_sys_opts), slave_yout)))

            result = (slave_tout, slave_yout, slave_flux, estimate_id[0], estimate_id[1], estimate_id[2],
                      perturbation_id)

        return data['task'], result


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
                                 'figure_format': 'eps',
                                 'ident_index_label': ['sample_name', 'data_set_id']})

        # retrieve identifiability data and process it for validation
        v1_ident.validation_info()

        user_ode_opts = {'iter': 'Newton', 'discr': 'Adams', 'atol': 1e-10, 'rtol': 1e-10,
                         'time_points': 200, 'display_progress': True, 'verbosity': 30}
        # initial ss to begin all simulations from
        y0 = np.array([5, 1, 1])
        # get and set true parameter values, if available separately
        default_parameters = true_parameter_values()

        v1_valid_obj = ValidateSim(kotte_model.kotte_ck_ode, kotte_model.kotte_ck_flux,
                                   **{'kinetics': 2, 'ode_opts': user_ode_opts, 't_final': 200, 'wt_y0': y0,
                                      'i_parameter': default_parameters, 'sample_size': 1, 'noise_std': 0.05,
                                      'validate_index_label': ['estimate_id', 'sample_name', 'data_set_id',
                                                               'experiment_id'],
                                      'validate_file_name': os.path.join(os.getcwd(), 'validate/v1_kcat_validate'),
                                      'original_exp_file': v1_ident.original_exp_file})

        parameter_estimates, estimate_info = v1_valid_obj.create_parameter_list(v1_ident.select_values)

        job = ParallelValidate(slaves=range(1, size))

        validate_results = job.run_all(task='initial_sim', **{'parameters': parameter_estimates,
                                                              'estimate_info': estimate_info, 'sim_obj': v1_valid_obj})
        job.terminate_slaves()

        # separate initial sims data from perturbation sims data and create validation df from dict
        v1_valid_obj.create_df(validate_results)

        # process validation info for plotting
        import pdb;pdb.set_trace()
        v1_valid_obj.process_validate()

        import pdb;pdb.set_trace()

    else:

        print('I am %s Slave with rank %s of %s' % (name, str(rank), str(size)))
        ValidateSlave().run()

    import pdb;pdb.set_trace()
    return None


if __name__ == '__main__':
    v1_validate()
    import pdb;pdb.set_trace()
    print('Done\n')


# def setup_parallel_validate(validate_sim_objects, perturbations):
#
#     name = MPI.Get_processor_name()
#     rank = MPI.COMM_WORLD.Get_rank()
#     size = MPI.COMM_WORLD.Get_size()
#
#     sim_result = {}
#     print('I am  %s rank %d (total %d)' % (name, rank, size))
#     if rank == 0:  # Master
#         validate_job = ParallelValidate(slaves=range(1, size))
#         # import pdb;pdb.set_trace()
#         sim_result = validate_job.run_i_value(validate_sim_objects, perturbations)
#
#         validate_job.terminate_slaves()
#     else:  # Any slave
#         MySlave().run()
#
#     # rearrange results based on order of parameters/initial values
#     # import pdb;pdb.set_trace()
#     return sim_result
