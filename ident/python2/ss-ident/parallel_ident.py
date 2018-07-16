from mpi4py import MPI
from mpi_master_slave import Master, Slave
from mpi_master_slave import WorkQueue
from identifiability_analysis import run_flux_ident
from run_ident import ModelIdent
from names_strings import true_parameter_values
import pandas as pd
import os.path
import kotte_model


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

        self.work_queue.add_work(data)

    def run_all(self, task, **kwargs):
        """parallel application core"""

        if task == 'ident':
            original_df = kwargs['exp_df']

            # remove experiment_id as index
            reset_df = original_df.reset_index('experiment_id')
            ident_fun = kwargs['ident_fun']
            flux_id = kwargs['flux_id']
            flux_choice = kwargs['flux_choice']

            idx = pd.IndexSlice
            all_df_indices = reset_df.index.unique().tolist()
            # create tuple of indices
            # import pdb; pdb.set_trace()
            for j_index, sample_data_set_id in enumerate(all_df_indices):
                j_exp_data_set = reset_df.loc[idx[sample_data_set_id],
                                              ['acetate', 'pep', 'fdp', 'E', 'v1', 'v2', 'v3', 'v5']].values.tolist()
                flat_data_list = [i_element for i_data in j_exp_data_set for i_element in i_data]
                self.__add_next_task(task=task, **{'exp_data': flat_data_list, 'ident_fun': ident_fun,
                                                   'flux_id': flux_id, 'flux_choice': flux_choice, 'sample_id':
                                                       sample_data_set_id[0], 'data_set_id': sample_data_set_id[1]})

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

        return results


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


def ident_processing(ident_object, ident_result):
    """call all required functions to process identifiability information"""
    # collect, arrange and collate data
    ordered_info = ident_object.order_ident_data(ident_result)
    ident_object.create_dict_for_df(ordered_info)

    # write all ident data to file
    ident_object.write_ident_info_file()

    # process ident info for further analysis
    ident_object.process_ident()

    # extract parameter va;ues for model validation
    ident_object.get_parameter_value()

    import pdb;pdb.set_trace()
    default_parameter_values = true_parameter_values()
    import pdb;pdb.set_trace()
    ident_object.parameter_values_plot(default_parameter_values, violin=True, box=False, bins=1)

    ident_object.identifiability_plot()

    ident_object.exp_info_plot()

    return None


def v1_kcat_ident():
    name = MPI.Get_processor_name()
    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()

    if rank == 0:  # master
        print('I am %s Master with rank %s of %s' % (name, str(rank), str(size)))
        v1_obj = ModelIdent(ident_fun=kotte_model.flux_1_kcat_ident,
                            arranged_data_file_name=os.path.join(os.getcwd(), 'exp/exp_v1_2_experiments'),
                            ident_data_file_name=os.path.join(os.getcwd(), 'ident/ident_v1_kcat'),
                            **{'original_exp_file': os.path.join(os.getcwd(), 'exp/experiments'),
                               'flux_id': 1, 'flux_choice': 2,
                               'values_figure': os.path.join(os.getcwd(), 'results/v1_kcat_parameter_values.eps'),
                               'ident_figure': os.path.join(os.getcwd(), 'results/v1_kcat_ident.eps'),
                               'exp_figure': os.path.join(os.getcwd(), 'results/v1_kcat_exp.eps'),
                               'figure_format': 'eps'})
        exp_df = v1_obj.retrieve_df_from_file()

        job = ParallelProcess(slaves=range(1, size))

        ident_result = job.run_all(task='ident', **{'exp_df': exp_df, 'ident_fun': v1_obj.ident_fun,
                                                    'flux_id': v1_obj.flux_id, 'flux_choice': v1_obj.flux_choice})

        job.terminate_slaves()

        # process ident data
        ident_processing(v1_obj, ident_result)

    else:
        print('I am %s Slave with rank %s of %s' % (name, str(rank), str(size)))
        ProcessSlave().run()

    return None


def v1_vmax_ident():
    name = MPI.Get_processor_name()
    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()

    if rank == 0:  # master
        print('I am %s Master with rank %s of %s' % (name, str(rank), str(size)))
        v1_obj = ModelIdent(ident_fun=kotte_model.flux_1_Vmax_ident,
                            arranged_data_file_name=os.path.join(os.getcwd(), 'exp/exp_v1_2_experiments'),
                            ident_data_file_name=os.path.join(os.getcwd(), 'ident/ident_v1_vmax'),
                            **{'original_exp_file': os.path.join(os.getcwd(), 'exp/experiments'),
                               'flux_id': 1, 'flux_choice': 1,
                               'values_figure': os.path.join(os.getcwd(), 'results/v1_vmax_parameter_values.eps'),
                               'ident_figure': os.path.join(os.getcwd(), 'results/v1_vmax_ident.eps'),
                               'exp_figure': os.path.join(os.getcwd(), 'results/v1_vmax_exp.eps'),
                               'figure_format': 'eps'})
        exp_df = v1_obj.retrieve_df_from_file()

        job = ParallelProcess(slaves=range(1, size))

        ident_result = job.run_all(task='ident', **{'exp_df': exp_df, 'ident_fun': v1_obj.ident_fun,
                                                    'flux_id': v1_obj.flux_id, 'flux_choice': v1_obj.flux_choice})

        job.terminate_slaves()

        # process ident data
        import pdb;pdb.set_trace()
        ident_processing(v1_obj, ident_result)

    else:
        print('I am %s Slave with rank %s of %s' % (name, str(rank), str(size)))
        ProcessSlave().run()

    return None


def v2_ident():
    name = MPI.Get_processor_name()
    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()

    if rank == 0:  # master
        print('I am %s Master with rank %s of %s' % (name, str(rank), str(size)))
        v2_obj = ModelIdent(ident_fun=kotte_model.flux_2_ident_expression,
                            arranged_data_file_name=os.path.join(os.getcwd(), 'exp/exp_v2_2_experiments'),
                            ident_data_file_name=os.path.join(os.getcwd(), 'ident/ident_v2'),
                            **{'original_exp_file': os.path.join(os.getcwd(), 'exp/experiments'),
                               'flux_id': 2, 'flux_choice': 0,
                               'values_figure': os.path.join(os.getcwd(), 'results/v2_parameter_values.eps'),
                               'ident_figure': os.path.join(os.getcwd(), 'results/v2_ident.eps'),
                               'exp_figure': os.path.join(os.getcwd(), 'results/v2_exp.eps'),
                               'figure_format': 'eps'})
        exp_df = v2_obj.retrieve_df_from_file()

        job = ParallelProcess(slaves=range(1, size))

        ident_result = job.run_all(task='ident', **{'exp_df': exp_df, 'ident_fun': v2_obj.ident_fun,
                                                    'flux_id': v2_obj.flux_id, 'flux_choice': v2_obj.flux_choice})

        job.terminate_slaves()

        # process ident data
        ident_processing(v2_obj, ident_result)

    else:
        print('I am %s Slave with rank %s of %s' % (name, str(rank), str(size)))
        ProcessSlave().run()

    return None


def v3_ident():
    name = MPI.Get_processor_name()
    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()

    if rank == 0:  # master
        print('I am %s Master with rank %s of %s' % (name, str(rank), str(size)))
        v3_obj = ModelIdent(ident_fun=kotte_model.flux_3_value1_ident,
                            arranged_data_file_name=os.path.join(os.getcwd(), 'exp/exp_v3_3_experiments'),
                            ident_data_file_name=os.path.join(os.getcwd(), 'ident/ident_v3'),
                            **{'original_exp_file': os.path.join(os.getcwd(), 'exp/experiments'),
                               'flux_id': 3, 'flux_choice': 1,
                               'values_figure': os.path.join(os.getcwd(), 'results/v3_parameter_values.eps'),
                               'ident_figure': os.path.join(os.getcwd(), 'results/v3_ident.eps'),
                               'exp_figure': os.path.join(os.getcwd(), 'results/v3_exp.eps'),
                               'figure_format': 'eps'})
        exp_df = v3_obj.retrieve_df_from_file()

        job = ParallelProcess(slaves=range(1, size))

        ident_result = job.run_all(task='ident', **{'exp_df': exp_df, 'ident_fun': v3_obj.ident_fun,
                                                    'flux_id': v3_obj.flux_id, 'flux_choice': v3_obj.flux_choice})

        job.terminate_slaves()

        # process ident data
        ident_processing(v3_obj, ident_result)

    else:
        print('I am %s Slave with rank %s of %s' % (name, str(rank), str(size)))
        ProcessSlave().run()

    return None


def v5_ident():
    name = MPI.Get_processor_name()
    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()

    if rank == 0:  # master
        print('I am %s Master with rank %s of %s' % (name, str(rank), str(size)))
        v5_obj = ModelIdent(ident_fun=kotte_model.flux_5_value2_ident,
                            arranged_data_file_name=os.path.join(os.getcwd(), 'exp/exp_v5_2_experiments'),
                            ident_data_file_name=os.path.join(os.getcwd(), 'ident/ident_v5'),
                            **{'original_exp_file': os.path.join(os.getcwd(), 'exp/experiments'),
                               'flux_id': 5, 'flux_choice': 2,
                               'values_figure': os.path.join(os.getcwd(), 'results/v5_parameter_values.eps'),
                               'ident_figure': os.path.join(os.getcwd(), 'results/v5_ident.eps'),
                               'exp_figure': os.path.join(os.getcwd(), 'results/v5_exp.eps'),
                               'figure_format': 'eps'})
        exp_df = v5_obj.retrieve_df_from_file()

        job = ParallelProcess(slaves=range(1, size))

        ident_result = job.run_all(task='ident', **{'exp_df': exp_df, 'ident_fun': v5_obj.ident_fun,
                                                    'flux_id': v5_obj.flux_id, 'flux_choice': v5_obj.flux_choice})

        job.terminate_slaves()

        # process ident data
        ident_processing(v5_obj, ident_result)

    else:
        print('I am %s Slave with rank %s of %s' % (name, str(rank), str(size)))
        ProcessSlave().run()

    return None


def v3_var_1_ident():
    name = MPI.Get_processor_name()
    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()

    if rank == 0:  # master
        print('I am %s Master with rank %s of %s' % (name, str(rank), str(size)))
        v3_obj = ModelIdent(ident_fun=kotte_model.flux_3_var1,
                            arranged_data_file_name=os.path.join(os.getcwd(), 'exp/exp_v3_2_experiments'),
                            ident_data_file_name=os.path.join(os.getcwd(), 'ident/ident_v3_k3fdp'),
                            **{'original_exp_file': os.path.join(os.getcwd(), 'exp/experiments'),
                               'flux_id': 6, 'flux_choice': 1,
                               'values_figure': os.path.join(os.getcwd(), 'results/v3_k3fdp_parameter_values.eps'),
                               'ident_figure': os.path.join(os.getcwd(), 'results/v3_k3fdp_ident.eps'),
                               'exp_figure': os.path.join(os.getcwd(), 'results/v3_k3fdp_exp.eps'),
                               'figure_format': 'eps'})
        exp_df = v3_obj.retrieve_df_from_file()

        job = ParallelProcess(slaves=range(1, size))

        ident_result = job.run_all(task='ident', **{'exp_df': exp_df, 'ident_fun': v3_obj.ident_fun,
                                                    'flux_id': v3_obj.flux_id, 'flux_choice': v3_obj.flux_choice})

        job.terminate_slaves()

        # process ident data
        ident_processing(v3_obj, ident_result)

    else:
        print('I am %s Slave with rank %s of %s' % (name, str(rank), str(size)))
        ProcessSlave().run()

    return None


def v3_var_2_ident():
    name = MPI.Get_processor_name()
    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()

    if rank == 0:  # master
        print('I am %s Master with rank %s of %s' % (name, str(rank), str(size)))
        v3_obj = ModelIdent(ident_fun=kotte_model.flux_3_var2,
                            arranged_data_file_name=os.path.join(os.getcwd(), 'exp/exp_v3_2_experiments'),
                            ident_data_file_name=os.path.join(os.getcwd(), 'ident/ident_v3_k3pep'),
                            **{'original_exp_file': os.path.join(os.getcwd(), 'exp/experiments'),
                               'flux_id': 6, 'flux_choice': 2,
                               'values_figure': os.path.join(os.getcwd(), 'results/v3_k3pep_parameter_values.eps'),
                               'ident_figure': os.path.join(os.getcwd(), 'results/v3_k3pep_ident.eps'),
                               'exp_figure': os.path.join(os.getcwd(), 'results/v3_k3pep_exp.eps'),
                               'figure_format': 'eps'})
        exp_df = v3_obj.retrieve_df_from_file()

        job = ParallelProcess(slaves=range(1, size))

        ident_result = job.run_all(task='ident', **{'exp_df': exp_df, 'ident_fun': v3_obj.ident_fun,
                                                    'flux_id': v3_obj.flux_id, 'flux_choice': v3_obj.flux_choice})

        job.terminate_slaves()

        # process ident data
        ident_processing(v3_obj, ident_result)

    else:
        print('I am %s Slave with rank %s of %s' % (name, str(rank), str(size)))
        ProcessSlave().run()

    return None


def v1_kcat_mwc_ident():
    name = MPI.Get_processor_name()
    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()

    if rank == 0:  # master
        print('I am %s Master with rank %s of %s' % (name, str(rank), str(size)))
        v1_obj = ModelIdent(ident_fun=kotte_model.flux_1_kcat_ident,
                            arranged_data_file_name=os.path.join(os.getcwd(), 'exp/exp_v1_2_experiments_mwc'),
                            ident_data_file_name=os.path.join(os.getcwd(), 'ident/ident_v1_kcat_mwc'),
                            **{'original_exp_file': os.path.join(os.getcwd(), 'exp/experiments_mwc'),
                               'flux_id': 1, 'flux_choice': 2,
                               'values_figure': os.path.join(os.getcwd(), 'results/v1_kcat_mwc_parameter_values.eps'),
                               'ident_figure': os.path.join(os.getcwd(), 'results/v1_kcat_mwc_ident.eps'),
                               'exp_figure': os.path.join(os.getcwd(), 'results/v1_kcat_mwc_exp.eps'),
                               'figure_format': 'eps'})
        exp_df = v1_obj.retrieve_df_from_file()

        job = ParallelProcess(slaves=range(1, size))

        ident_result = job.run_all(task='ident', **{'exp_df': exp_df, 'ident_fun': v1_obj.ident_fun,
                                                    'flux_id': v1_obj.flux_id, 'flux_choice': v1_obj.flux_choice})

        job.terminate_slaves()

        # process ident data
        ident_processing(v1_obj, ident_result)

    else:
        print('I am %s Slave with rank %s of %s' % (name, str(rank), str(size)))
        ProcessSlave().run()

    return None


def v1_vmax_mwc_ident():
    name = MPI.Get_processor_name()
    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()

    if rank == 0:  # master
        print('I am %s Master with rank %s of %s' % (name, str(rank), str(size)))
        v1_obj = ModelIdent(ident_fun=kotte_model.flux_1_Vmax_ident,
                            arranged_data_file_name=os.path.join(os.getcwd(), 'exp/exp_v1_2_experiments_mwc'),
                            ident_data_file_name=os.path.join(os.getcwd(), 'ident/ident_v1_vmax_mwc'),
                            **{'original_exp_file': os.path.join(os.getcwd(), 'exp/experiments_mwc'),
                               'flux_id': 1, 'flux_choice': 1,
                               'values_figure': os.path.join(os.getcwd(), 'results/v1_vmax_mwc_parameter_values.eps'),
                               'ident_figure': os.path.join(os.getcwd(), 'results/v1_vmax_mwc_ident.eps'),
                               'exp_figure': os.path.join(os.getcwd(), 'results/v1_vmax_mwc_exp.eps'),
                               'figure_format': 'eps'})
        exp_df = v1_obj.retrieve_df_from_file()

        job = ParallelProcess(slaves=range(1, size))

        ident_result = job.run_all(task='ident', **{'exp_df': exp_df, 'ident_fun': v1_obj.ident_fun,
                                                    'flux_id': v1_obj.flux_id, 'flux_choice': v1_obj.flux_choice})

        job.terminate_slaves()

        # process ident data
        ident_processing(v1_obj, ident_result)

    else:
        print('I am %s Slave with rank %s of %s' % (name, str(rank), str(size)))
        ProcessSlave().run()
    return None


def v2_mwc_ident():
    name = MPI.Get_processor_name()
    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()

    if rank == 0:  # master
        print('I am %s Master with rank %s of %s' % (name, str(rank), str(size)))
        v2_obj = ModelIdent(ident_fun=kotte_model.flux_2_ident_expression,
                            arranged_data_file_name=os.path.join(os.getcwd(), 'exp/exp_v2_2_experiments_mwc'),
                            ident_data_file_name=os.path.join(os.getcwd(), 'ident/ident_v2_mwc'),
                            **{'original_exp_file': os.path.join(os.getcwd(), 'exp/experiments_mwc'),
                               'flux_id': 2, 'flux_choice': 0,
                               'values_figure': os.path.join(os.getcwd(), 'results/v2_mwc_parameter_values.eps'),
                               'ident_figure': os.path.join(os.getcwd(), 'results/v2_mwc_ident.eps'),
                               'exp_figure': os.path.join(os.getcwd(), 'results/v2_mwc_exp.eps'),
                               'figure_format': 'eps'})
        exp_df = v2_obj.retrieve_df_from_file()

        job = ParallelProcess(slaves=range(1, size))

        ident_result = job.run_all(task='ident', **{'exp_df': exp_df, 'ident_fun': v2_obj.ident_fun,
                                                    'flux_id': v2_obj.flux_id, 'flux_choice': v2_obj.flux_choice})

        job.terminate_slaves()

        # process ident data
        ident_processing(v2_obj, ident_result)

    else:
        print('I am %s Slave with rank %s of %s' % (name, str(rank), str(size)))
        ProcessSlave().run()
    return None


def v3_mwc_ident():
    name = MPI.Get_processor_name()
    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()

    if rank == 0:  # master
        print('I am %s Master with rank %s of %s' % (name, str(rank), str(size)))
        v3_obj = ModelIdent(ident_fun=kotte_model.flux_3_value1_ident,
                            arranged_data_file_name=os.path.join(os.getcwd(), 'exp/exp_v3_3_experiments_mwc'),
                            ident_data_file_name=os.path.join(os.getcwd(), 'ident/ident_v3_mwc'),
                            **{'original_exp_file': os.path.join(os.getcwd(), 'exp/experiments_mwc'),
                               'flux_id': 3, 'flux_choice': 1,
                               'values_figure': os.path.join(os.getcwd(), 'results/v3_mwc_parameter_values.eps'),
                               'ident_figure': os.path.join(os.getcwd(), 'results/v3_mwc_ident.eps'),
                               'exp_figure': os.path.join(os.getcwd(), 'results/v3_mwc_exp.eps'),
                               'figure_format': 'eps'})
        exp_df = v3_obj.retrieve_df_from_file()

        job = ParallelProcess(slaves=range(1, size))

        ident_result = job.run_all(task='ident', **{'exp_df': exp_df, 'ident_fun': v3_obj.ident_fun,
                                                    'flux_id': v3_obj.flux_id, 'flux_choice': v3_obj.flux_choice})

        job.terminate_slaves()

        # process ident data
        ident_processing(v3_obj, ident_result)

    else:
        print('I am %s Slave with rank %s of %s' % (name, str(rank), str(size)))
        ProcessSlave().run()
    return None


def v5_mwc_ident():
    name = MPI.Get_processor_name()
    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()

    if rank == 0:  # master
        print('I am %s Master with rank %s of %s' % (name, str(rank), str(size)))
        v5_obj = ModelIdent(ident_fun=kotte_model.flux_5_value2_ident,
                            arranged_data_file_name=os.path.join(os.getcwd(), 'exp/exp_v5_2_experiments_mwc'),
                            ident_data_file_name=os.path.join(os.getcwd(), 'ident/ident_v5_mwc'),
                            **{'original_exp_file': os.path.join(os.getcwd(), 'exp/experiments_mwc'),
                               'flux_id': 5, 'flux_choice': 2,
                               'values_figure': os.path.join(os.getcwd(), 'results/v5_mwc_parameter_values.eps'),
                               'ident_figure': os.path.join(os.getcwd(), 'results/v5_mwc_ident.eps'),
                               'exp_figure': os.path.join(os.getcwd(), 'results/v5_mwc_exp.eps'),
                               'figure_format': 'eps'})
        exp_df = v5_obj.retrieve_df_from_file()

        job = ParallelProcess(slaves=range(1, size))

        ident_result = job.run_all(task='ident', **{'exp_df': exp_df, 'ident_fun': v5_obj.ident_fun,
                                                    'flux_id': v5_obj.flux_id, 'flux_choice': v5_obj.flux_choice})

        job.terminate_slaves()

        # process ident data
        ident_processing(v5_obj, ident_result)

    else:
        print('I am %s Slave with rank %s of %s' % (name, str(rank), str(size)))
        ProcessSlave().run()
    return None


if __name__ == '__main__':
    v1_kcat_ident()
    v1_vmax_ident()
    import pdb;pdb.set_trace()
    v2_ident()
    v3_ident()
    v3_var_1_ident()
    v3_var_2_ident()
    v5_ident()
    v1_kcat_mwc_ident()
    v1_vmax_mwc_ident()
    v2_mwc_ident()
    v3_mwc_ident()
    v5_mwc_ident()
    import pdb;pdb.set_trace()
    print('Done\n')