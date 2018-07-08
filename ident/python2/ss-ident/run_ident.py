from parallel_ident import setup_parallel_ident
from names_strings import default_ident_parameter_name, true_parameter_values
from collections import defaultdict
from plot_ident_results import plot_on_axis_object_box, plot_on_axis_object_hist, plot_on_axis_object_violin
from plot_ident_results import plot_on_axis_object, set_hbar_axis_properties
from plot_ident_results import plot_on_axis_object_polar, set_polar_axis_limits
import pandas as pd
import itertools as it
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.font_manager as fnt
import matplotlib as mpl
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

    def get_parameter_value(self):
        """extract parameter values in a given flux to re-simulate model with newly determined parameters.
        get parameter values from data sets that can detect all parameters"""

        import pdb;pdb.set_trace()
        # retrieve ident data from file
        ident_df = self.retrieve_ident_df_from_file()

        # get data sets identifying each parameter
        identifying_data_sets = [set(i_parameter_data_set) for i_parameter_data_set in self.processed_info["sample_data_set_id"]]
        size_of_data_sets = [len(i_parameter_set) for i_parameter_set in identifying_data_sets]
        sort_index = np.argsort(size_of_data_sets)  # last position is the biggest data set
        largest_set = identifying_data_sets[sort_index[-1]]
        # del identifying_data_sets[sort_index[-1]]
        for i_index in range(len(sort_index) - 2, -1, -1):
            largest_set.intersection_update(identifying_data_sets[sort_index[i_index]])

        # get parameter values from data sets in largest_set
        idx = pd.IndexSlice
        relevant_df = ident_df.loc[idx[list(largest_set)], ['parameter_name', 'parameter_value']]
        all_parameter_values = []
        all_parameter_names = []
        all_parameter_data_sets = [list(largest_set) for _ in self.processed_info["parameter_names"]]
        for i_parameter, i_parameter_name in enumerate(self.processed_info["parameter_names"]):
            i_p_value = relevant_df[relevant_df["parameter_name"] == i_parameter_name]["parameter_value"].values
            all_parameter_values.append(i_p_value)
            all_parameter_names.append(i_parameter_name)
        all_parameter_info = {"names": all_parameter_names, "values": all_parameter_values,
                              "data_sets": all_parameter_data_sets,
                              "total_data_sets": self.processed_info["total_data_sets"],
                              "flux_name": self.processed_info["flux_name"]}

        import pdb;pdb.set_trace()
        return all_parameter_info

    def parameter_values_plot(self, original_values=(), violin=False, box=True, bins=1):
        """plot distribution of parameter values as a box plot, violin plot and/or histogram"""
        number_parameters = len(self.processed_info["parameter_names"])
        if box:
            f1 = plt.figure(figsize=(10, 8), dpi=100)
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
            box_axis.set_ylabel('Parameter values')
            box_axis.tick_params(axis='both', which='major', direction='in', length=3.5, width=0.5, color='black',
                                 bottom=True)

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
                hist_axis.set_xlabel('Parameter value')
                hist_axis.tick_params(axis='both', which='major', direction='in', length=3.5, width=0.5, color='black',
                                      bottom=True)
            # f1.savefig(self.ident_figure, dpi=300, format='pdf', facecolor='w', edgecolor='k', transparent=True)
            f1.savefig(self.values_figure, format=self.figure_format, transparent=True, frameon=True,
                       bbox_inches='tight')

        if violin:
            f2 = plt.figure(figsize=(10, 8), dpi=100)
            plot_grid = gridspec.GridSpec(2, number_parameters)

            # plot box plot
            violin_axis = f2.add_subplot(plot_grid[0, :])
            plot_on_axis_object_violin(violin_axis, self.processed_info["parameter_values"])
            violin_axis.set_xticks(np.arange(1, len(self.processed_info['parameter_names']) + 1))
            violin_axis.set_xticklabels(self.processed_info['parameter_names'])
            violin_axis.set_ylabel('Parameter values')
            # violin_axis.grid(b=False)
            violin_axis.tick_params(axis='both', which='major', direction='in', length=3.5, width=0.5, color='black',
                                    bottom=True)

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
                hist_axis.set_xlabel('Parameter value')
                # hist_axis.grid(b=False)
                hist_axis.tick_params(axis='both', which='major', direction='in', length=3.5, width=0.5, color='black',
                                      bottom=True)

            f2.savefig(self.values_figure, format=self.figure_format, transparent=True, frameon=True,
                       bbox_inches='tight')
        return None

    def identifiability_plot(self):
        """plot identifiability (number and percentage of data sets identifying each parameter)
        of every parameter in a given flux"""
        number_of_parameters = len(self.processed_info["parameter_names"])
        f, ax = plt.subplots(1, 1, figsize=(10, 8), dpi=100)
        x_data = self.processed_info["ident_mean"]
        y_data = range(0, number_of_parameters)
        gap = 0.05
        # y_add = [0, gap + 0.2, 2 * (gap + 0.2)]
        y_add = [i_y_data * (gap + 0.2) for i_y_data in y_data]
        x_error = self.processed_info["ident_std"]
        x_percent_mean = self.processed_info["ident_percent_mean"]
        x_percent_std = self.processed_info["ident_percent_std"]
        x_max = max(self.processed_info["total_data_sets"])
        x_label = ['Number of data combinations used for identification']
        plot_on_axis_object(ax, x_data, y_add, x_error=x_error)
        set_hbar_axis_properties(ax, y_add, y_tick_label=self.processed_info["parameter_names"], x_max=x_max,
                                 x_percent_mean=x_percent_mean, x_percent_std=x_percent_std, x_label=x_label,
                                 figure_title=self.processed_info["flux_name"][0] + ' parameters')
        ax.set_axis_on()
        ax.tick_params(axis='both', which='major', direction='in', length=3.5, width=0.5, color='black', bottom=True)
        f.savefig(self.ident_figure, format=self.figure_format, transparent=False, frameon=True,
                  bbox_inches='tight')
        return None

    def exp_info_plot(self):
        """plot experiment contribution frequency for each parameter in a polar plot"""
        mpl.rcParams['lines.linewidth'] = 2
        mpl.rcParams['lines.color'] = 'r'
        number_parameters = len(self.processed_info["parameter_names"])
        if number_parameters >= 3:
            number_of_columns = 3
        else:
            number_of_columns = 1
        if number_parameters % number_of_columns != 0:
            number_of_rows = (number_parameters + 1) / number_of_columns
        else:
            number_of_rows = number_parameters / number_of_columns
        # figure = plt.figure(figsize=(6, 4))
        # inner_grid = gridspec.GridSpec(number_of_rows, number_of_columns, wspace=0.2, hspace=0.2)
        f, ax = plt.subplots(int(number_of_rows), int(number_of_columns), subplot_kw=dict(projection='polar'),
                             figsize=(10, 8), dpi=100, gridspec_kw={"wspace": 0.2, "hspace": 0.2})
        all_max_y_data = []
        for i_parameter, i_parameter_info in enumerate(self.processed_info["exp_info"]):
            # ax = plt.Subplot(figure, inner_grid[i_parameter])
            fill_colors = ['b', 'g', 'y', 'r']
            # plot data from all positions
            pos_labels = []
            for i_pos, (i_pos_key, i_pos_val) in enumerate(i_parameter_info.items()):
                pos_labels.append(i_pos_key)
                y_data = i_pos_val["frequency"]
                x_labels = i_pos_val["names"]
                max_y_data = plot_on_axis_object_polar(ax[i_parameter], x_data=x_labels, y_data=y_data,
                                                       data_label=pos_labels[-1], fill_color=fill_colors[i_pos])
                ax[i_parameter].set_axis_on()
                all_max_y_data.append(max_y_data)

        # set axis limits for all polar plots on subplot
        for i_parameter, (_, i_parameter_name) in enumerate(zip(self.processed_info["exp_info"],
                                                                self.processed_info["parameter_names"])):
            set_polar_axis_limits(ax[i_parameter], max(all_max_y_data))
            f_p = fnt.FontProperties(size=14, weight='demibold')
            ax[i_parameter].set_title(i_parameter_name, **{'font_properties': f_p})

        # add legend
        plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))

        f.savefig(self.exp_figure, format=self.figure_format, transparent=True, frameon=True,
                  bbox_inches='tight')
        return None


if __name__ == '__main__':
    # extract experimental data from file
    default_parameter_values = true_parameter_values()

    v1_ident = ModelIdent(ident_fun=kotte_model.flux_1_kcat_ident,
                          arranged_data_file_name=os.path.join(os.getcwd(), 'exp/exp_v1_2_experiments'),
                          ident_data_file_name=os.path.join(os.getcwd(), 'ident/ident_v1_kcat'),
                          **{'original_exp_file': os.path.join(os.getcwd(), 'exp/experiments'),
                             'flux_id': 1, 'flux_choice': 2,
                             'values_figure': os.path.join(os.getcwd(), 'results/v1_kcat_parameter_values.eps'),
                             'ident_figure': os.path.join(os.getcwd(), 'results/v1_kcat_ident.eps'),
                             'exp_figure': os.path.join(os.getcwd(), 'results/v1_kcat_exp.eps'),
                             'figure_format': 'eps'})
    # test identifiability
    print('Practical Identifiability Analysis of v1 with 2 parameters: k1cat and K1ac\n')
    ident_data_df = v1_ident.perform_ident()

    v1_ident.process_ident()

    v1_ident.get_parameter_value()

    # get parameter value plot
    v1_ident.parameter_values_plot(default_parameter_values, violin=True, box=False, bins=1)

    v1_ident.identifiability_plot()

    v1_ident.exp_info_plot()

    import pdb;pdb.set_trace()
    print('Done\n')

    # # test identifiability and store data to file
    # from kotte_model import flux_ident_2_data_combination
    # flux_ident_2_data_combination(arranged_data_df, flux_ids=[1], flux_choice=[2],
    #                               ident_fun_choice=ident_fun_choice, file_name=storage_file_name)