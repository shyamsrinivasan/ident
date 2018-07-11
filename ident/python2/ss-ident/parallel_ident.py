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


def v1_ident():
    name = MPI.Get_processor_name()
    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()

    if rank == 0:  # master
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

        # collect, arrange and collate data
        ordered_info = v1_obj.order_ident_data(ident_result)
        v1_obj.create_dict_for_df(ordered_info)

        # write all ident data to file
        v1_obj.write_ident_info_file()

        # process ident info for further analysis
        v1_obj.process_ident()

        # extract parameter va;ues for model validation
        v1_obj.get_parameter_value()

        default_parameter_values = true_parameter_values()
        v1_obj.parameter_values_plot(default_parameter_values, violin=True, box=False, bins=1)

        v1_obj.identifiability_plot()

        v1_obj.exp_info_plot()
    else:
        ProcessSlave().run()

    return v1_obj


def v2_ident():
    name = MPI.Get_processor_name()
    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()

    if rank == 0:  # master
        v2_obj = ModelIdent(ident_fun=kotte_model.flux_2_ident_expression,
                            arranged_data_file_name=os.path.join(os.getcwd(), 'exp/exp_v2_2_experiments'),
                            ident_data_file_name=os.path.join(os.getcwd(), 'ident/ident_v2'),
                            **{'original_exp_file': os.path.join(os.getcwd(), 'exp/experiments'),
                               'flux_id': 1, 'flux_choice': 1,
                               'values_figure': os.path.join(os.getcwd(), 'results/v2_parameter_values.eps'),
                               'ident_figure': os.path.join(os.getcwd(), 'results/v2_ident.eps'),
                               'exp_figure': os.path.join(os.getcwd(), 'results/v2_exp.eps'),
                               'figure_format': 'eps'})
        exp_df = v2_obj.retrieve_df_from_file()

        job = ParallelProcess(slaves=range(1, size))

        ident_result = job.run_all(task='ident', **{'exp_df': exp_df, 'ident_fun': v2_obj.ident_fun,
                                                    'flux_id': v2_obj.flux_id, 'flux_choice': v2_obj.flux_choice})

        job.terminate_slaves()

        # collect, arrange and collate data
        ordered_info = v2_obj.order_ident_data(ident_result)
        v2_obj.create_dict_for_df(ordered_info)

        # write all ident data to file
        v2_obj.write_ident_info_file()

        # process ident info for further analysis
        v2_obj.process_ident()

        # extract parameter va;ues for model validation
        v2_obj.get_parameter_value()

        default_parameter_values = true_parameter_values()
        v2_obj.parameter_values_plot(default_parameter_values, violin=True, box=False, bins=1)

        v2_obj.identifiability_plot()

        v2_obj.exp_info_plot()
    else:
        ProcessSlave().run()

    return v2_obj


if __name__ == '__main__':

    # v1_object = v1_ident()
    
    v2_object = v2_ident()
    import pdb;pdb.set_trace()
    print('Done\n')