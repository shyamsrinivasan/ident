from mpi4py import MPI
from mpi_master_slave import Master, Slave
from mpi_master_slave import WorkQueue
from identifiability_analysis import run_flux_ident
from simulate_ode import simulate_ode
from names_strings import default_ident_parameter_name, true_parameter_values
import pandas as pd
import numpy as np
import kotte_model
import os.path


class ModelIdent(object):
    def __init__(self, ident_fun, arranged_data_file_name, ident_data_file_name, **kwargs):
        self.ident_fun = ident_fun
        self.arranged_data_file = arranged_data_file_name
        self.ident_file = ident_data_file_name

        try:
            self.original_exp_file = kwargs['original_exp_file']
        except KeyError:
            self.original_exp_file = ''
        try:
            self.flux_id = kwargs['flux_id']
            self.flux_name = 'flux{}'.format(self.flux_id)
        except KeyError:
            self.flux_id = []
            self.flux_name = ''
        try:
            self.flux_choice = kwargs['flux_choice']
        except KeyError:
            self.flux_choice = []
        try:
            self.parameter_name = default_ident_parameter_name(flux_name=self.flux_name, flux_choice=self.flux_choice)
        except AttributeError:
            self.parameter_name = []

        try:
            self.arranged_index_label = kwargs['arranged_index_label']
        except KeyError:
            self.arranged_index_label = ['sample_name', 'data_set_id', 'experiment_id']
        try:
            self.original_index_label = kwargs['original_index_label']
        except KeyError:
            self.original_index_label = ['sample_name', 'experiment_id']

        try:
            self.figure_format = kwargs['figure_format']
        except KeyError:
            self.figure_format = 'eps'

        try:
            self.values_figure = kwargs['values_figure']
        except KeyError:
            self.values_figure = ''

        try:
            self.ident_figure = kwargs['ident_figure']
        except KeyError:
            self.ident_figure = ''

        try:
            self.exp_figure = kwargs['exp_figure']
        except KeyError:
            self.exp_figure = ''

        self.unique_indices = []
        self.ident_data = {}
        self.ident_index_label = []
        self.processed_info = {}
        self.select_values = {}

    def retrieve_df_from_file(self, original_exp=0):
        """retrieve experimental data from csv file"""
        # read dataframe from csv file
        experiment_df = pd.read_csv(self.arranged_data_file, index_col=self.arranged_index_label)

        # lexographic ordering of exp df indices
        experiment_df.sort_index(level='sample_name', inplace=True)
        experiment_df.sort_index(level='data_set_id', inplace=True)

        if original_exp:
            experiment_df = pd.read_csv(self.original_exp_file, index_col=self.original_index_label)

        return experiment_df

    def retrieve_ident_df_from_file(self):
        """retrieve identifiability data from csv file to dataframe"""
        # read dataframe from csv file
        ident_df = pd.read_csv(self.ident_file, index_col=self.ident_index_label)

        # lexicographic sorting of indices in ident_df
        ident_df.sort_index(level='sample_name', inplace=True)
        ident_df.sort_index(level='data_set_id', inplace=True)

        return ident_df


class ParallelProcess(object):
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
        if task == 'ident':
            data = {'ident_fun': kwargs['ident_fun'], 'flux_id': kwargs['flux_id'],
                    'flux_choice': kwargs['flux_choice'], 'sample_id': kwargs['sample_id'],
                    'data_set_id': kwargs['data_set_id'], 'exp_data': kwargs['exp_data'], 'task': task}
            # data = [ident_fun, ident_args]

        elif task == 'initial_sim':
            data = {'ode_fun': kwargs['ode_rhs_fun'], 'y0': kwargs['initial_value'], 'id': kwargs['estimate_id'],
                    'ode_sys_opts': kwargs['ode_sys_opts'], 'ode_opts': kwargs['ode_opts'],
                    't_final': kwargs['t_final'], 'flux_fun': kwargs['flux_fun'], 'task': task}

        elif task == 'validate':
            # set ode rhs function, initial condition and parameters
            data = {'ode_fun': kwargs['ode_rhs_fun'], 'y0': kwargs['initial_value'], 'id': kwargs['estimate_id'],
                    'ode_sys_opts': kwargs['ode_sys_opts'], 'ode_opts': kwargs['ode_opts'],
                    't_final': kwargs['t_final'], 'flux_fun': kwargs['flux_fun'],
                    'perturbations': kwargs['perturbations'], 'task': task}
            # data = [validate_fun, validate_args]
        else:
            data = {}

        self.work_queue.add_work(data)

    def run_all(self, task, **kwargs):
        """parallel application core"""
        if task == 'ident':
            exp_df = kwargs['exp_df']
            ident_fun = kwargs['ident_fun']
            flux_id = kwargs['flux_id']
            flux_choice = kwargs['flux_choice']

            idx = pd.IndexSlice
            all_df_indices = exp_df.index.unique().tolist()
            # create tuple of indices
            # import pdb; pdb.set_trace()
            for j_index, sample_data_set_id in enumerate(all_df_indices):
                j_exp_data_set = exp_df.loc[idx[sample_data_set_id],
                                            ['acetate', 'pep', 'fdp', 'E', 'v1', 'v2', 'v3', 'v5']].values.tolist()
                flat_data_list = [i_element for i_data in j_exp_data_set for i_element in i_data]
                self.__add_next_task(task=task, **{'exp_data': flat_data_list, 'ident_fun': ident_fun,
                                                   'flux_id': flux_id, 'flux_choice': flux_choice, 'sample_id':
                                                       sample_data_set_id[0], 'data_set_id': sample_data_set_id[1]})

        elif task == 'initial_sim':
            estimates = kwargs['parameters']
            estimate_id = kwargs['estimate_info']
            sim_obj = kwargs['sim_obj']
            for j_index, (j_estimate, j_estimate_id) in enumerate(zip(estimates, estimate_id)):
                self.__add_next_task(task=task, **{'ode_rhs_fun': sim_obj.ode_rhs_fun,
                                                   'flux_fun': sim_obj.flux_fun, 't_final': sim_obj.t_final,
                                                   'parameters': j_estimate, 'ode_opts': sim_obj.ode_opts,
                                                   'initial_value': sim_obj.y0, 'estimate_id': j_estimate_id})

        elif task == 'validate_sim':
            pass

        # Keeep starting slaves as long as there is work to do
        results = []
        while not self.work_queue.done():
            # give more work to do to each idle slave (if any)
            self.work_queue.do_work()

            # reclaim returned data from completed slaves
            for slave_return_data in self.work_queue.get_completed_work():
                task, output = slave_return_data
                if task == 'ident':
                    ident_info, slave_flux_id, slave_flux_choice, sample_id, data_set_id = output
                    i_slave_result = {'flux_id': slave_flux_id, 'flux_choice': slave_flux_choice, 'sample_id': sample_id,
                                      'data_set_id': data_set_id, 'ident_info': ident_info}
                    results.append(i_slave_result)

                    print('Master: slave finished task %s returning: %s)' % (task, str(data_set_id)))

                elif task == 'initial_sim':
                    time_course, y_result, flux, estimate_id, sample_id, data_set_id = output
                    i_slave_result = {'t': time_course, 'y': y_result, 'flux': flux, 'estimate_id': estimate_id,
                                      'sample_id': sample_id, 'data_set_id': data_set_id}
                    results.append(i_slave_result)

        return results


    # def run_ident(self, exp_df, ident_fun, flux_id, flux_choice):
    #     """
    #     This is the core where I keep starting slaves
    #     as long as there is work to do
    #     """
    #
    #     # let's prepare our work queue. This can be built at initialization time
    #     # but it can also be added later as more work become available
    #     #
    #     idx = pd.IndexSlice
    #     all_df_indices = exp_df.index.unique().tolist()
    #     # create tuple of indices
    #     # import pdb; pdb.set_trace()
    #     for j_index, sample_data_set_id in enumerate(all_df_indices):
    #         j_exp_data_set = exp_df.loc[idx[sample_data_set_id],
    #                                     ['acetate', 'pep', 'fdp', 'E', 'v1', 'v2', 'v3', 'v5']].values.tolist()
    #         flat_data_list = [i_element for i_data in j_exp_data_set for i_element in i_data]
    #         self.__add_next_task(task='ident', **{'exp_data': flat_data_list, 'ident_fun': ident_fun,
    #                                               'flux_id': flux_id, 'flux_choice': flux_choice, 'sample_id':
    #                                                   sample_data_set_id[0], 'data_set_id': sample_data_set_id[1]})
    #
    #     # Keeep starting slaves as long as there is work to do
    #     ident_results = []
    #     while not self.work_queue.done():
    #         # give more work to do to each idle slave (if any)
    #         self.work_queue.do_work()
    #         # reclaim returned data from completed slaves
    #         for slave_return_data in self.work_queue.get_completed_work():
    #             task, output = slave_return_data
    #             if task == 'ident':
    #                 ident_info, slave_flux_id, slave_flux_choice, sample_id, data_set_id = output
    #                 i_slave_result = {'flux_id': flux_id, 'flux_choice': flux_choice, 'sample_id': sample_id,
    #                                   'data_set_id': data_set_id, 'ident_info': ident_info}
    #                 ident_results.append(i_slave_result)
    #
    #                 print('Master: slave finished its task returning: %s)' % str(data_set_id))
    #         # sleep some time
    #         # time.sleep(0.01)
    #     return ident_results

    # def run_initial_sim(self, sim_obj, estimates, estimate_id):
    #     """run initial simulation from given y0 for all given parameter sets"""
    #     for j_index, (j_estimate, j_estimate_id) in enumerate(zip(estimates, estimate_id)):
    #         self.__add_next_task(task='initial_sim', **{'ode_rhs_fun': sim_obj.ode_rhs_fun,
    #                                                     'flux_fun': sim_obj.flux_fun, 't_final': sim_obj.t_final,
    #                                                     'parameters': j_estimate, 'ode_opts': sim_obj.ode_opts,
    #                                                     'initial_value': sim_obj.y0, 'estimate_id': j_estimate_id})
    #
    #     # Keeep starting slaves as long as there is work to do
    #     all_boolean = []
    #     all_tout = []
    #     all_yout = []
    #     all_y0_id = []
    #     while not self.work_queue.done():
    #         # give more work to do to each idle slave (if any)
    #         self.work_queue.do_work()
    #         # reclaim returned data from completed slaves
    #         for slave_return_data in self.work_queue.get_completed_work():
    #             y0_id, done, time_course, y_result = slave_return_data
    #             # import pdb; pdb.set_trace()
    #             all_boolean.append(done)
    #             all_tout.append(time_course)
    #             all_yout.append(y_result)
    #             all_y0_id.append(y0_id)
    #
    #             if done:
    #                 print('Master: slave finished its task returning: %s)' % str(y0_id))
    #         # sleep some time
    #         # time.sleep(0.3)
    #     results = {'time': all_tout, 'y': all_yout, 'id': all_y0_id, 'boolean': all_boolean}
    #     return results
    #     # Keeep starting slaves as long as there is work to do
    #     ident_results = []
    #     while not self.work_queue.done():
    #         # give more work to do to each idle slave (if any)
    #         self.work_queue.do_work()
    #         # reclaim returned data from completed slaves
    #         for slave_return_data in self.work_queue.get_completed_work():
    #             task, output = slave_return_data
    #             if task == 'ident':
    #                 ident_info, slave_flux_id, slave_flux_choice, sample_id, data_set_id = output
    #                 i_slave_result = {'flux_id': flux_id, 'flux_choice': flux_choice, 'sample_id': sample_id,
    #                                   'data_set_id': data_set_id, 'ident_info': ident_info}
    #                 ident_results.append(i_slave_result)
    #
    #                 print('Master: slave finished its task returning: %s)' % str(data_set_id))
    #         # sleep some time
    #         # time.sleep(0.01)
    #     return ident_results


    # def run_validation(self, sim_object, estimates, estimate_ids):
    #     """
    #     This is the core where I keep starting slaves
    #     as long as there is work to do
    #     """
    #     # let's prepare our work queue. This can be built at initialization time
    #     # but it can also be added later as more work become available
    #
    #     # import pdb; pdb.set_trace()
    #     for j_index, sample_data_set_id in enumerate(all_df_indices):
    #         j_exp_data_set = ident_df.loc[idx[sample_data_set_id],
    #                                     ['acetate', 'pep', 'fdp', 'E', 'v1', 'v2', 'v3', 'v5']].values.tolist()
    #         flat_data_list = [i_element for i_data in j_exp_data_set for i_element in i_data]
    #         self.__add_next_task(task='validate', **{'exp_data': flat_data_list, 'ident_fun': ident_fun,
    #                                               'flux_id': flux_id, 'flux_choice': flux_choice, 'sample_id':
    #                                                   sample_data_set_id[0], 'data_set_id': sample_data_set_id[1]})
    #
    #     # Keeep starting slaves as long as there is work to do
    #     ident_results = []
    #     while not self.work_queue.done():
    #         # give more work to do to each idle slave (if any)
    #         self.work_queue.do_work()
    #         # reclaim returned data from completed slaves
    #         for slave_return_data in self.work_queue.get_completed_work():
    #             task, output = slave_return_data
    #             if task == 'validate':
    #                 ident_info, slave_flux_id, slave_flux_choice, sample_id, data_set_id = output
    #                 i_slave_result = {'flux_id': flux_id, 'flux_choice': flux_choice, 'sample_id': sample_id,
    #                                   'data_set_id': data_set_id, 'ident_info': ident_info}
    #                 ident_results.append(i_slave_result)
    #
    #                 print('Master: slave finished its task returning: %s)' % str(data_set_id))
    #         # sleep some time
    #         # time.sleep(0.01)
    #     return validate_results


class ProcessSlave(Slave):
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
        super(ProcessSlave, self).__init__()

    def do_work(self, data):
        """do work method overrides Slave.do_work() and defines work
        to be done by every slave"""
        rank = MPI.COMM_WORLD.Get_rank()
        name = MPI.Get_processor_name()

        print('  Slave %s rank %d executing task %s' % (name, rank, data['task']))
        if data['task'] == 'ident':
            ident_fun = data['ident_fun']
            exp_data = data['exp_data']
            flux_id = data['flux_id']
            flux_choice = data['flux_choice']
            sample_id = data['sample_id']
            data_set_id = data['data_set_id']

            print('  Slave %s rank %d executing ident for sample: %s, data set: %s' %
                  (name, rank, sample_id, data_set_id))

            ident_info, _, _ = run_flux_ident(ident_fun, exp_data, flux_id, flux_choice)

            return data['task'], (ident_info, flux_id, flux_choice, sample_id, data_set_id)

        elif data['task'] == 'initial_sim':
            # define explicit assimulo problem
            rhs_fun = data['ode_fun']
            y_initial = data['y0']
            estimate_id = data['estimate_id']
            ode_opts = data['ode_opts']
            ode_sys_opts = data['ode_sys_opts']
            t_final = data['t_final']
            all_options = [ode_opts, ode_sys_opts]

            print('  Slave %s rank %d executing ident for estimate: %s sample: %s, data set: %s' %
                  (name, rank, estimate_id[0], estimate_id[1], estimate_id[2]))
            time_course, y_result, _, _ = simulate_ode(rhs_fun, y_initial, tf=t_final, opts=all_options)

            # calculate flux
            flux_fun = data['flux_fun']
            flux = np.array(list(map(lambda x: flux_fun(x, ode_sys_opts), y_result)))

            return data['task'], (time_course, y_result, flux, estimate_id[0], estimate_id[1], estimate_id[2])

        elif data['task'] == 'validate':
            pass

            return data['task'], ([], [], [], [], [])


# def setup_parallel_ident(ident_fun, flux_id, flux_choice, exp_data, name, rank, size):
#     # name = MPI.Get_processor_name()
#     # rank = MPI.COMM_WORLD.Get_rank()
#     # size = MPI.COMM_WORLD.Get_size()
#
#     sim_result = {}
#     print('I am  %s rank %d (total %d)' % (name, rank, size))
#     if rank == 0:  # Master
#         ident_job = ParallelProcess(slaves=range(1, size))
#         # import pdb;pdb.set_trace()
#         sim_result = ident_job.run_all(exp_data, ident_fun=ident_fun, flux_id=flux_id, flux_choice=flux_choice)
#
#         ident_job.terminate_slaves()
#     else:  # Any slave
#         MySlave().run()
#
#     return sim_result


if __name__ == '__main__':
    # user_ode_opts = {'iter': 'Newton', 'discr': 'Adams', 'atol': 1e-10, 'rtol': 1e-10,
    #                  'time_points': 200, 'display_progress': True, 'verbosity': 30}
    # # initial ss to begin all simulations from
    # y0 = np.array([5, 1, 1])
    # # get and set true parameter values, if available separately
    # default_parameters = true_parameter_values()
    # # create simulation object to simulate model with above parameters and initial conditions
    # model_1 = ModelSim(kotte_model.kotte_ck_ode, kotte_model.kotte_ck_flux, noise=0, **{'kinetics': 2,
    #                                                                                     'ode_opts': user_ode_opts,
    #                                                                                     't_final': 200,
    #                                                                                     'wt_y0': y0,
    #                                                                                     'i_parameter':
    #                                                                                         default_parameters})
    name = MPI.Get_processor_name()
    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()

    if rank == 0:  # master
        v1_ident = ModelIdent(ident_fun=kotte_model.flux_1_kcat_ident,
                              arranged_data_file_name=os.path.join(os.getcwd(), 'exp/exp_v1_2_experiments'),
                              ident_data_file_name=os.path.join(os.getcwd(), 'ident/ident_v1_kcat'),
                              **{'original_exp_file': os.path.join(os.getcwd(), 'exp/experiments'),
                                 'flux_id': 1, 'flux_choice': 2,
                                 'values_figure': os.path.join(os.getcwd(), 'results/v1_kcat_parameter_values.eps'),
                                 'ident_figure': os.path.join(os.getcwd(), 'results/v1_kcat_ident.eps'),
                                 'exp_figure': os.path.join(os.getcwd(), 'results/v1_kcat_exp.eps'),
                                 'figure_format': 'eps'})
        exp_df = v1_ident.retrieve_df_from_file()

        job = ParallelProcess(slaves=range(1, size))
        ident_result = job.run_all(task='ident', **{'exp_df': exp_df, 'ident_fun': v1_ident.ident_fun,
                                                    'flux_id': v1_ident.flux_id, 'flux_choice': v1_ident.flux_choice})
    else:
        ProcessSlave.run()
