from parallel_ident import setup_parallel_ident
from names_strings import default_ident_parameter_name
from collections import defaultdict
import pandas as pd
import itertools as it
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

        self.unique_indices = []
        self.ident_data = []
        self.ident_index_label = []

    def retrieve_df_from_file(self, original_exp=0):
        """retrieve experimental data from csv file"""
        # read dataframe from csv file
        experiment_df = pd.read_csv(self.arranged_data_file, index_col=self.arranged_index_label)

        if original_exp:
            experiment_df = pd.read_csv(self.original_exp_file, index_col=self.original_index_label)
        return experiment_df

    def retrieve_ident_df_from_file(self):
        """retrieve identifiability data from csv file to dataframe"""
        # read dataframe from csv file
        ident_df = pd.read_csv(self.ident_file, index_col=self.ident_index_label)
        return ident_df

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
        ordered_info = self.__order_ident_data(sim_result)
        self.__create_dict_for_df(ordered_info)
        final_result_df = self.__write_ident_info_file()

        return final_result_df

    def __order_ident_data(self, all_results):
        """collect ident info for all data sets from all samples together"""

        # read experimental data from file (serially)
        arranged_df = self.retrieve_df_from_file()

        reset_df = arranged_df.reset_index('experiment_id')
        # lexographic ordering of df indices
        reset_df.sort_index(level='data_set_id', inplace=True)
        reset_df.sort_index(level='sample_name', inplace=True)

        # arrange results based on original experimental df
        all_df_indices = reset_df.index.unique().tolist()
        self.unique_indices = all_df_indices
        ordered_info = [j_result for i_data_id in self.unique_indices for j_result in all_results
                        if (j_result['sample_id'] == i_data_id[0] and j_result['data_set_id'] == i_data_id[1])]

        return ordered_info

    def __create_dict_for_df(self, ident_results):
        """create dictionary of identifiability results for future writing to file"""
        # read experimental data from file (serially)
        arranged_df = self.retrieve_df_from_file()
        # remove index experiment_id
        reset_df = arranged_df.reset_index('experiment_id')
        # lexographic ordering of remaining df indices
        reset_df.sort_index(level='data_set_id', inplace=True)
        reset_df.sort_index(level='sample_name', inplace=True)
        temp_dict = {}
        all_data = defaultdict(list)
        empty_dict = {}

        for i_data_set in ident_results:
            data_set_id = i_data_set['data_set_id']
            sample_id = i_data_set['sample_id']
            for i_parameter, i_parameter_info in enumerate(i_data_set['ident_info']):
                temp_dict["flux_name"] = self.flux_name
                temp_dict["flux_choice"] = self.flux_choice
                temp_dict["parameter_name"] = self.parameter_name[i_parameter]
                i_parameter_nr, i_parameter_dr, i_parameter_value = i_parameter_info
                temp_dict["parameter_nr"] = i_parameter_nr
                temp_dict["parameter_dr"] = i_parameter_dr
                temp_dict["parameter_value"] = i_parameter_value
                temp_dict["data_set_id"] = data_set_id
                temp_dict["sample_name"] = sample_id
                if i_parameter_value > 0:
                    temp_dict["identified"] = True
                else:
                    temp_dict["identified"] = False
                for key, value in it.chain(empty_dict.items(), temp_dict.items()):
                    all_data[key].append(value)
        self.ident_data = all_data
        return self

    def __write_ident_info_file(self):
        """create data frame from identifiability data and write to csv file for future use"""
        # read experimental data from file (serially)
        arranged_df = self.retrieve_df_from_file()

        # reset index of experimental data df
        reset_exp_df = arranged_df.reset_index("sample_name")
        reset_exp_df.sort_index(level='data_set_id', inplace=True)
        reset_exp_df.reset_index("data_set_id", inplace=True)

        # create data frame
        ident_df = pd.DataFrame(self.ident_data, columns=self.ident_data.keys())

        # number of occurrences of each data set id = number of experiments per data set in the first sample
        number_samples = len(ident_df["sample_name"].unique())
        first_sample_rows = ident_df[ident_df["sample_name"] == 'sample_0']
        data_set_id_frequency = int(max(first_sample_rows["data_set_id"].value_counts()))

        # all experiment ids
        experiment_pos_names = ['experiment_{}_id'.format(i_experiment) for i_experiment in
                                range(0, data_set_id_frequency)]
        experiment_pos_parameters = ['experiment_{}_parameter'.format(i_experiment)
                                     for i_experiment in range(0, data_set_id_frequency)]

        # extract experiment ids for each data set
        # get all data set ids
        data_set_ids = reset_exp_df["data_set_id"].unique()

        # get experiments for each data set based on first sample only
        all_data_set_experiments = [reset_exp_df[(reset_exp_df["sample_name"] == "sample_0") &
                                                 (reset_exp_df["data_set_id"] == j_data_set_id)]
                                    ["parameter_name"].index.values.tolist() for j_data_set_id in data_set_ids]
        all_data_set_exp_parameters = [reset_exp_df[(reset_exp_df["sample_name"] == "sample_0") &
                                                    (reset_exp_df["data_set_id"] == j_data_set_id)]
                                       ["parameter_name"].values.tolist() for j_data_set_id in data_set_ids]
        all_pos_experiment_id = [[i_p for j_data_set in all_data_set_experiments
                                  for i_p in [j_data_set[j_position_exp]] * len(experiment_pos_names)] * number_samples
                                 for j_position_exp in range(0, len(experiment_pos_names))]
        all_pos_exp_parameters = [[i_p for j_data_set in all_data_set_exp_parameters
                                   for i_p in
                                   [j_data_set[j_position_exp]] * len(experiment_pos_parameters)] * number_samples
                                  for j_position_exp in range(0, len(experiment_pos_parameters))]
        experiment_pos_info_keys = experiment_pos_names + experiment_pos_parameters
        experiment_pos_info_values = all_pos_experiment_id + all_pos_exp_parameters
        exp_info_dict = dict(zip(experiment_pos_info_keys, experiment_pos_info_values))

        self.ident_data.update(exp_info_dict)

        # multi index tuples
        ind_tuple = [(j_sample, j_data_set) for j_sample, j_data_set in
                     zip(self.ident_data["sample_name"], self.ident_data["data_set_id"])]

        # multi index index
        self.ident_index_label = ['sample_name', 'data_set_id']
        index = pd.MultiIndex.from_tuples(ind_tuple, names=self.ident_index_label)

        # remove redundant columns
        temp_ident_data = self.ident_data
        del temp_ident_data["sample_name"]
        del temp_ident_data["data_set_id"]

        # create multi index data frame
        ident_data_df = pd.DataFrame(temp_ident_data, index=index, columns=temp_ident_data.keys())

        # save data frame to csv file
        ident_data_df.to_csv(self.ident_file, index_label=self.ident_index_label)
        print('Identifiability Data written to given file\n')
        return ident_data_df


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
    ident_data_df = v1_ident.perform_ident()

    # # test identifiability and store data to file
    # from kotte_model import flux_ident_2_data_combination
    # flux_ident_2_data_combination(arranged_data_df, flux_ids=[1], flux_choice=[2],
    #                               ident_fun_choice=ident_fun_choice, file_name=storage_file_name)