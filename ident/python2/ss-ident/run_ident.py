from parallel_ident import setup_parallel_ident
from names_strings import default_ident_parameter_name, true_parameter_values
from collections import defaultdict
from plot_ident_results import plot_on_axis_object_box, plot_on_axis_object_hist, plot_on_axis_object_violin
import pandas as pd
import itertools as it
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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

        try:
            self.ident_figure = kwargs['ident_figure']
        except KeyError:
            self.ident_figure = ''

        self.unique_indices = []
        self.ident_data = {}
        self.ident_index_label = []
        self.processed_info = {}

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

    def perform_ident(self):
        # read experimental data from file (serially)
        arranged_df = self.retrieve_df_from_file()
        reset_df = arranged_df.reset_index('experiment_id')

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

        # arrange results based on original experimental df
        all_df_indices = reset_df.index.unique().tolist()
        self.unique_indices = all_df_indices
        ordered_info = [j_result for i_data_id in self.unique_indices for j_result in all_results
                        if (j_result['sample_id'] == i_data_id[0] and j_result['data_set_id'] == i_data_id[1])]

        return ordered_info

    def __create_dict_for_df(self, ident_results):
        """create dictionary of identifiability results for future writing to file"""
        # read experimental data from file (serially)
        # arranged_df = self.retrieve_df_from_file()
        # remove index experiment_id
        # reset_df = arranged_df.reset_index('experiment_id')
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

    @staticmethod
    def __parameter_ident_info(ident_df):
        """returns info on identifiability and value of each parameter in df"""
        # get parameter names
        all_parameter_names = ident_df["parameter_name"].unique().tolist()

        all_p_values = []
        all_p_data_set_names = []
        all_p_identifiability = []
        total_data_sets = []
        all_p_sample_names = []
        for i_parameter_name in all_parameter_names:
            # get all data sets identifying each parameter
            identifying_df = ident_df[(ident_df["parameter_name"] == i_parameter_name) & (ident_df["identified"])]

            # get data set names
            all_p_data_set_names.append([i_value[1] for i_value in identifying_df.index.values])
            all_p_identifiability.append(len(all_p_data_set_names[-1]))

            # get parameter values
            all_p_values.append([np.array(i_value) for i_value in identifying_df["parameter_value"].values])
            total_data_sets.append(len(ident_df.index.levels[1]))
            all_p_sample_names.append([i_value[0] for i_value in identifying_df.index.values])

        all_p_info = {"parameter_names": all_parameter_names,
                      "parameter_values": all_p_values,
                      "identifiability": all_p_identifiability,
                      "data_sets": all_p_data_set_names,
                      "total_data_sets": total_data_sets,
                      "identifiability_percentage": [np.array(float(i_nr) * 100 / float(i_dr)) for i_nr, i_dr in
                                                     zip(all_p_identifiability, total_data_sets)],
                      "sample_name": all_p_sample_names}
        return all_p_info

    @staticmethod
    def __sample_ident_info(all_sample_info):
        """return calculated mean identifiability and std in identifiability
        from all samples for all parameters"""
        # get common data sets identifying all parameters in all samples
        number_parameters = len(all_sample_info[0]["parameter_names"])

        all_identifiabilites = [i_sample_info["identifiability"] for i_sample_info in all_sample_info]
        sample_mean_identifiability = np.mean(np.array(all_identifiabilites), axis=0)
        sample_std_identifiability = np.std(np.array(all_identifiabilites), axis=0)

        sample_mean_ident = [np.array(i_parameter_ident) for i_parameter_ident in sample_mean_identifiability]
        sample_std_ident = [np.array(i_parameter_ident) for i_parameter_ident in sample_std_identifiability]

        all_identifiabilites_percent = [i_sample_info["identifiability_percentage"]
                                        for i_sample_info in all_sample_info]

        sample_mean_identifiability_percent = np.mean(np.array(all_identifiabilites_percent), axis=0)
        sample_std_identifiability_percent = np.std(np.array(all_identifiabilites_percent), axis=0)

        sample_mean_ident_percent = [np.array(i_parameter_ident) for i_parameter_ident in
                                     sample_mean_identifiability_percent]
        sample_std_ident_percent = [np.array(i_parameter_ident) for i_parameter_ident in
                                    sample_std_identifiability_percent]

        all_sample_data_pair = [[i_p for i_sample in all_sample_info
                                 for i_p in
                                 zip(i_sample["sample_name"][i_parameter], i_sample["data_sets"][i_parameter])]
                                for i_parameter in range(0, number_parameters)]
        ident_dict = {"ident_mean": sample_mean_ident,
                      "ident_std": sample_std_ident,
                      "ident_percent_mean": sample_mean_ident_percent,
                      "ident_percent_std": sample_std_ident_percent,
                      "sample_data_set_id": all_sample_data_pair}
        return ident_dict

    @staticmethod
    def __parameter_exp_info(ident_df, exp_df, parameter_ident_info):
        """return experiment type information for each parameter"""
        # get parameter names
        all_parameter_names = ident_df["parameter_name"].unique().tolist()
        number_experiments = len(all_parameter_names)
        exp_column_ids = ['experiment_{}_parameter'.format(i_experiment) for i_experiment in
                          range(0, number_experiments)]
        exp_column_name = ['experiment_{}'.format(i_experiment) for i_experiment in range(0, number_experiments)]

        # all possible ss perturbation experiment types classified on the basis of parameter perturbed
        all_possible_perturbations = set.union(set(exp_df["parameter_name"].unique()),
                                               {'wt', 'ac', 'k1cat', 'V2max', 'V3max'})
        all_parameter_exp_info = []
        for i_parameter, i_parameter_name in enumerate(all_parameter_names):
            # get all data sets identifying each parameter
            identifying_df = ident_df[(ident_df["parameter_name"] == i_parameter_name) & (ident_df["identified"])]
            all_experiment_info = {}

            # get frequency of each experiment
            for i_experiment, i_experiment_pos in enumerate(exp_column_ids):
                exp_frequency = identifying_df[i_experiment_pos].value_counts()

                # get name value pairs
                name_value_pair = [(j_name, np.array(float(j_value) * 100 / parameter_ident_info[i_parameter]))
                                   for j_name, j_value in zip(exp_frequency.index.values, exp_frequency.values)]

                # add missing experiment type with value = 0
                missing_perturbation = all_possible_perturbations.difference(exp_frequency.index.values)
                name_value_pair = name_value_pair + list(zip(missing_perturbation, [0.0] * len(missing_perturbation)))

                # arrange name/values in desired order for every parameter
                given_value = []
                for i_given_name in all_possible_perturbations:
                    given_value.append([i_obtained_value for i_obtained_name, i_obtained_value in name_value_pair
                                        if i_given_name == i_obtained_name][0])
                all_experiment_info.update({exp_column_name[i_experiment]: {"names": list(all_possible_perturbations),
                                                                            "frequency": given_value}})
            all_parameter_exp_info.append(all_experiment_info)

        return all_parameter_exp_info

    def process_ident(self, ident=1, exp_info=1):
        """process ident data to create final data frame of results for plotting"""

        # retrieve arranged experimental data from file
        exp_df = self.retrieve_df_from_file()

        # lexographic ordering of exp df indices
        exp_df.sort_index(level='experiment_id', inplace=True)

        # retrieve ident data from file
        ident_df = self.retrieve_ident_df_from_file()

        idx = pd.IndexSlice

        # number of samples
        sample_names = ident_df.index.levels[0].values.tolist()
        number_samples = len(sample_names)

        if ident:
            # parameter identifiability information for each sample
            all_sample_info = []
            for i_sample in sample_names:
                j_sample_info = ident_df.loc[idx[i_sample, :], :]
                # get all data sets identifying each parameter in each sample
                all_p_info = self.__parameter_ident_info(j_sample_info)
                all_sample_info.append(all_p_info)

            ident_dict = self.__sample_ident_info(all_sample_info)
        else:
            ident_dict = {}

        # parameter identifiability information for all samples along with mean and std between samples
        all_parameter_info = self.__parameter_ident_info(ident_df)
        all_parameter_info.update(ident_dict)

        # get experiment information for each parameter (only for noise-less data)
        # number_experiments = len(all_parameter_info["names"])
        if exp_info and number_samples == 1:
            exp_info = self.__parameter_exp_info(ident_df, exp_df, all_parameter_info["ident_mean"])
        else:
            exp_info = []
        all_parameter_info.update({"exp_info": exp_info})

        # get flux names
        all_flux_names = ident_df["flux_name"].unique().tolist() * len(all_parameter_info["parameter_names"])
        all_parameter_info.update({"flux_name": all_flux_names})

        self.processed_info = all_parameter_info

        return self

    def parameter_values_plot(self, original_values=(), violin=False, box=True, bins=1):
        """plot distribution of parameter values as a box plot, violin plot and/or histogram"""
        number_parameters = len(self.processed_info["parameter_names"])
        if box:
            f1 = plt.figure(figsize=(10, 8), dpi=100, tight_layout=True)
            plot_grid = gridspec.GridSpec(2, number_parameters)

            # plot box plot
            box_axis = f1.add_subplot(plot_grid[0, :])
            if original_values:
                plot_on_axis_object_box(box_axis, self.processed_info["parameter_values"],
                                        mark_value=[original_values[i_name] for i_name in
                                                    self.processed_info["parameter_names"]])
            else:
                plot_on_axis_object_box(box_axis, self.processed_info["parameter_values"],
                                        mark_value=[])
            box_axis.set_xticklabels(self.processed_info['parameter_names'])

            # plot histogram
            for i_parameter, (i_parameter_value, i_parameter_name) in enumerate(
                    zip(self.processed_info["parameter_values"], self.processed_info['parameter_names'])):
                # parameter_name = info_dict["names"][i_parameter]
                hist_axis = f1.add_subplot(plot_grid[1, i_parameter])
                if original_values:
                    plot_on_axis_object_hist(hist_axis, i_parameter_value, mark_value=original_values[i_parameter_name],
                                             parameter_name=i_parameter_name, bins=bins)
                else:
                    plot_on_axis_object_hist(hist_axis, i_parameter_value, mark_value=[],
                                             parameter_name=i_parameter_name, bins=bins)

            f1.savefig(self.ident_figure, dpi=1000, format='pdf')

        if violin:
            f2 = plt.figure(figsize=(10, 8), dpi=100, tight_layout=True)
            plot_grid = gridspec.GridSpec(2, number_parameters)

            # plot box plot
            violin_axis = f2.add_subplot(plot_grid[0, :])
            plot_on_axis_object_violin(violin_axis, self.processed_info["parameter_values"])
            violin_axis.set_xticks(np.arange(1, len(self.processed_info['parameter_names']) + 1))
            violin_axis.set_xticklabels(self.processed_info['parameter_names'])

            # plot histogram
            for i_parameter, (i_parameter_value, i_parameter_name) in enumerate(
                    zip(self.processed_info["parameter_values"], self.processed_info['parameter_names'])):
                # parameter_name = info_dict["names"][i_parameter]
                hist_axis = f2.add_subplot(plot_grid[1, i_parameter])
                if original_values:
                    plot_on_axis_object_hist(hist_axis, i_parameter_value, mark_value=original_values[i_parameter_name],
                                             parameter_name=i_parameter_name)
                else:
                    plot_on_axis_object_hist(hist_axis, i_parameter_value, mark_value=[],
                                             parameter_name=i_parameter_name)

            f2.savefig(self.ident_figure, dpi=1000, format='pdf')
        return None


if __name__ == '__main__':
    # extract experimental data from file
    exp_figure = os.path.join(os.getcwd(), 'results/v1_kcat_exp.eps')
    default_parameter_values = true_parameter_values()

    v1_ident = ModelIdent(ident_fun=kotte_model.flux_1_kcat_ident,
                          arranged_data_file_name=os.path.join(os.getcwd(), 'exp/exp_v1_2_experiments'),
                          ident_data_file_name=os.path.join(os.getcwd(), 'ident/ident_v1_kcat'),
                          **{'original_exp_file': os.path.join(os.getcwd(), 'exp/experiments'),
                             'flux_id': 1, 'flux_choice': 2,
                             'ident_figure': os.path.join(os.getcwd(), 'results/v1_kcat_ident')})
    # test identifiability
    print('Practical Identifiability Analysis of v1 with 2 parameters: k1cat and K1ac\n')
    ident_data_df = v1_ident.perform_ident()

    v1_ident.process_ident()

    # get parameter value plot
    v1_ident.parameter_values_plot(default_parameter_values, violin=True, box=False, bins=1)

    # # test identifiability and store data to file
    # from kotte_model import flux_ident_2_data_combination
    # flux_ident_2_data_combination(arranged_data_df, flux_ids=[1], flux_choice=[2],
    #                               ident_fun_choice=ident_fun_choice, file_name=storage_file_name)