import pandas as pd
from parallel_ident import setup_parallel_ident
import os.path
import kotte_model


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
        except KeyError:
            self.flux_id = []
        try:
            self.flux_choice = kwargs['flux_choice']
        except KeyError:
            self.flux_choice = []
        try:
            self.arranged_index_label = kwargs['arranged_index_label']
        except KeyError:
            self.arranged_index_label = ['sample_name', 'data_set_id', 'experiment_id']
        try:
            self.original_index_label = kwargs['original_index_label']
        except KeyError:
            self.original_index_label = ['sample_name', 'experiment_id']

    def retrieve_df_from_file(self, original_exp=0):
        """retrieve experimental data from csv file"""
        # read dataframe from csv file
        experiment_df = pd.read_csv(self.arranged_data_file, index_col=self.arranged_index_label)

        if original_exp:
            experiment_df = pd.read_csv(self.arranged_data_file, index_col=self.original_index_label)
        return experiment_df

    def arrange_ident_data(self):
        pass
        return self

    def perform_ident(self):
        # read experimental data from file (serially)
        arranged_df = self.retrieve_df_from_file()
        
        reset_df = arranged_df.reset_index('experiment_id')
        # lexographic ordering of df indices
        reset_df.sort_index(level='data_set_id', inplace=True)
        reset_df.sort_index(level='sample_name', inplace=True)

        # run ident analysis in parallel
        sim_result = setup_parallel_ident(ident_fun=self.ident_fun, flux_id=self.flux_id, flux_choice=self.flux_choice,
                                          exp_data=reset_df)
        # collect, arrange and collate data
        import pdb; pdb.set_trace()
        return sim_result

    # @staticmethod
    # def collect_data(exp_df, j_sample, numerical=0):
    #     """collect data from df from all data sets in single sample to
    #     peform identiiability analysis with cas and numerical method"""
    #     idx = pd.IndexSlice
    #     all_sample_ids = exp_df.index.levels[0].tolist()
    #     all_data_set_ids = exp_df.index.levels[1].tolist()
    #     # number_data_sets = (len(all_data_sets))
    #     all_exp_data = []
    #     if numerical:
    #         # avoid converting list of lists to single list
    #         for j_data_set, data_set_id in enumerate(all_data_set_ids):
    #             ident_data = exp_df.loc[idx[j_sample, data_set_id],
    #                                     ['acetate', 'pep', 'fdp', 'E', 'v1', 'v2', 'v3', 'v5']].values.tolist()
    #             single_list = [i_exp_data for i_exp_data in ident_data]
    #             all_exp_data.append(single_list)
    #     else:
    #         for j_data_set, data_set_id in enumerate(all_data_set_ids):
    #             ident_data = exp_df.loc[idx[j_sample, data_set_id],
    #                                     ['acetate', 'pep', 'fdp', 'E', 'v1', 'v2', 'v3', 'v5']].values.tolist()
    #             single_list = [i_variable for i_exp_data in ident_data for i_variable in i_exp_data]
    #             all_exp_data.append(single_list)
    #     return all_exp_data, all_data_set_ids

    # def multi_sample_ident_fun(self, ident_fun_list, all_data_df, flux_ids, flux_choice):
    #     """perform identifibaility analysis for multiple samples by
    #     looping through each experimental data sample for a list of identifibaility functions"""
    #
    #     reset_df = all_data_df.reset_index('experiment_id')
    #     # lexographic ordering of df indices
    #     reset_df.sort_index(level='data_set_id', inplace=True)
    #     reset_df.sort_index(level='sample_name', inplace=True)
    #
    #
    #     sample_ids = list(reset_df.index.levels[0])
    #     number_samples = (len(sample_ids))
    #     all_sample_ident_details = []
    #     for i_sample, i_sample_id in enumerate(sample_ids):
    #         print('Identifiability analysis with Data Sample Number {} of {}\n'.format(i_sample, number_samples))
    #         # collect experimental data from all data sets
    #         all_exp_data, _ = collect_data(reset_df, i_sample_id)
    #         # run identifiability with i_sample_data
    #         all_ident_values, _ = get_ident_value(ident_fun_list, all_exp_data, flux_ids, flux_choice)
    #         all_sample_ident_details.append(all_ident_values)
    #
    #     # initial information processing to get dictionary of relevant info for each flux and each parameter
    #     all_data = defaultdict(list)
    #     empty_dict = {}
    #     for i_sample, sample_data in enumerate(all_sample_ident_details):
    #         sample_name = 'sample_{}'.format(i_sample)
    #         # number_data_sets = len(sample_data)
    #         all_data = collect_ident_data(sample_name, sample_data, flux_ids, flux_choice, all_data, empty_dict)
    #
    #     return all_data


if __name__ == '__main__':
    # extract experimental data from file
    original_experiment_file = os.path.join(os.getcwd(), 'exp/experiments')
    arranged_data_file = os.path.join(os.getcwd(), 'exp/exp_v1_2_experiments')
    storage_file_name = os.path.join(os.getcwd(), 'ident/ident_v1_kcat')

    v1_ident = ModelIdent(ident_fun=kotte_model.flux_1_kcat_ident, arranged_data_file_name=arranged_data_file,
                          ident_data_file_name=storage_file_name, **{'original_exp_file': original_experiment_file,
                                                                     'flux_id': 1, 'flux_choice': 2})
    # test identifiability
    print('Practical Identifiability Analysis of v1 with 2 parameters: k1cat and K1ac\n')
    v1_ident.perform_ident()

    # # test identifiability and store data to file
    # from kotte_model import flux_ident_2_data_combination
    # flux_ident_2_data_combination(arranged_data_df, flux_ids=[1], flux_choice=[2],
    #                               ident_fun_choice=ident_fun_choice, file_name=storage_file_name)