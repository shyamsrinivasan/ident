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

        self.unique_indices = []

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
        return sim_result

    def collect_ident_data(self, all_results, all_data_dict, empty_dict):
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
        import pdb;pdb.set_trace()
        collated_info = [j_result for i_data_id in self.unique_indices for j_result in all_results
                         if (j_result['sample_id'] == i_data_id[0] and j_result['data_set_id'] == i_data_id[1])]

        # temp_dict = {}
        # for j_data_set, j_data_set_info in enumerate(j_sample_ident_data):
        #     data_set_name = 'data_set_{}'.format(j_data_set)
        #     for j_flux, j_flux_data in enumerate(j_data_set_info):
        #         flux_name = 'flux{}'.format(flux_ids[j_flux])
        #         # all_parameter names
        #         all_parameter_names = [ident_parameter_name(j_parameter,
        #                                                     flux_name,
        #                                                     flux_choice[j_flux])
        #                                for j_parameter in range(0, len(j_flux_data))]
        #         for i_parameter, i_parameter_info in enumerate(j_flux_data):
        #             temp_dict["flux_name"] = flux_name
        #             temp_dict["flux_choice"] = flux_choice[j_flux]
        #             # replace with call to parameter name file
        #             temp_dict["parameter_name"] = all_parameter_names[i_parameter]
        #             i_parameter_nr, i_parameter_dr, i_parameter_value = i_parameter_info
        #             temp_dict["parameter_nr"] = i_parameter_nr
        #             temp_dict["parameter_dr"] = i_parameter_dr
        #             temp_dict["parameter_value"] = i_parameter_value
        #             temp_dict["data_set_id"] = data_set_name
        #             temp_dict["sample_name"] = j_sample_name
        #             if i_parameter_value > 0:
        #                 temp_dict["identified"] = True
        #             else:
        #                 temp_dict["identified"] = False
        #             for key, value in it.chain(empty_dict.items(), temp_dict.items()):
        #                 all_data_dict[key].append(value)
        import pdb;pdb.set_trace()
        return collated_info





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