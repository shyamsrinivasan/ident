from run_sims import ModelSim
from names_strings import variable_name
from plot_ident_results import plot_on_axis_object_box, plot_on_axis_object_violin
from plot_ident_results import validation_hist, validation_scatter
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np


class ValidateSim(ModelSim):
    def __init__(self, rhs_fun, flux_fun, noise=0, **kwargs):
        super(ValidateSim, self).__init__(rhs_fun, flux_fun, noise, **kwargs)
        # self.rhs_fun = rhs_fun
        # self.flux_fun = flux_fun
        # self.noise = noise

        # all parameter perturbations to be tested
        try:
            self.test_perturbations = kwargs['test_perturbations']
        except KeyError:
            self.test_perturbations = [{"wt": 0}, {"ac": 1}, {"ac": 4}, {"ac": 9}, {"ac": -.1}, {"ac": -.5},
                                       {"k1cat": .1}, {"k1cat": .5}, {"k1cat": 1}, {"k1cat": -.1}, {"k1cat": -.5},
                                       {"V3max": .1}, {"V3max": .5}, {"V3max": 1}, {"V3max": -.1}, {"V3max": -.5},
                                       {"V2max": .1}, {"V2max": .5}, {"V2max": 1}, {"V2max": -.1}, {"V2max": -.5}]

        # self.wt_ss = []
        # self.wt_dynamic = []
        self.perturbation_ss = {}
        # self.perturbation_dynamic = []
        self.estimated_parameters = []
        self.estimate_ids = []

        try:
            self.validate_index_label = kwargs['validate_index_label']
        except KeyError:
            self.validate_index_label = []

        try:
            self.validate_file = kwargs['validate_file_name']
        except KeyError:
            self.validate_file = []

        try:
            self.original_exp_file = kwargs['original_exp_file']
        except KeyError:
            self.original_exp_file = []

        try:
            self.original_index_label = kwargs['original_exp_label']
        except KeyError:
            self.original_index_label = ['sample_name', 'experiment_id']

        self.concentration_validation = {}
        self.flux_validation = {}

    def get_name_value_parameter_pairs(self, estimate_info):
        """get estimated parameter values as name value pairs"""
        number_estimates = len(estimate_info['data_sets'])
        parameter_name_value_pair = [dict(zip(estimate_info['parameter_names'],
                                              [estimate_info['parameter_values'][i_parameter][i_estimate]
                                               for i_parameter, _ in enumerate(estimate_info['parameter_names'])]))
                                     for i_estimate in range(0, number_estimates)]
        self.estimated_parameters = parameter_name_value_pair
        return None

    def create_parameter_list(self, estimate_info):
        """create name value pairs of estimated parameters followed by list of all parameters for use in validation"""
        # create dictionary (of length n_p) of parameters
        self.get_name_value_parameter_pairs(estimate_info)

        # create list of all parameter values of size n_p with each of the above estimated values
        parameter_list = [self.i_parameter for _ in self.estimated_parameters]
        # data_set_id = []
        estimate_data_set_info = []
        for i_index, i_value in enumerate(parameter_list):
            # data_set_id.append(estimate_info['data_sets'][i_index])
            # estimate_id.append('estimate_{}'.format(i_index))
            estimate_data_set_info.append(('estimate_{}'.format(i_index), estimate_info['data_sets'][i_index][0],
                                           estimate_info['data_sets'][i_index][1]))
            for i_key in self.estimated_parameters[i_index].keys():
                i_value[i_key] = self.estimated_parameters[i_index][i_key]

        return parameter_list, estimate_data_set_info

    @staticmethod
    def separate_initial_perturbation(all_sim_results):
        """separate initial ss simulations from perturbation simulation
        for further processing of perturbation sims"""

        initial_sims = []
        perturbation_sims = []

        for j_result in all_sim_results:
            if j_result['initial']:
                initial_sims.append(j_result)
            elif j_result['perturbation']:
                perturbation_sims.append(j_result)

        return initial_sims, perturbation_sims

    @staticmethod
    def separate_ss_dyn(all_results):
        """get ss and dynamic information separately"""

        ss_info = []
        dynamic_info = []

        for j_result in all_results:
            ss_info.append(j_result['ss'])
            dynamic_info.append(j_result['dynamic'])

        return ss_info, dynamic_info

    @staticmethod
    def __create_ss_dict(ss_info, variable_type, noise=0):
        """create dictionary of ss values (concentration/fluxes)"""
        number_variables = len(ss_info[0])
        variable_name_info = [variable_name(variable_type, j_variable) for j_variable in range(0, number_variables)]

        variable_value_info = []

        for j_variable in range(0, number_variables):
            j_variable_info = []
            for i_experiment_id, i_experiment_info in enumerate(ss_info):
                j_variable_info.append(i_experiment_info[j_variable])
            variable_value_info.append(j_variable_info)

        return variable_name_info, variable_value_info

    def __convert_to_ss_dict_for_df(self, ss_info):
        """convert ss information to dictionary for creating df and writing to file"""

        # create dict of perturbation_ss values for writing to df and file
        # concentrations
        y_names, y_values = self.__create_ss_dict(ss_info=[i_ss['y'] for i_ss in ss_info],
                                                  variable_type='metabolite', noise=self.noise)
        y_dict = dict(zip(y_names, y_values))

        # fluxes
        f_names, f_values = self.__create_ss_dict(ss_info=[i_ss['flux'] for i_ss in ss_info],
                                                  variable_type='flux', noise=self.noise)
        f_dict = dict(zip(f_names, f_values))

        # get stable ss information
        final_ss_id = [j_ss_info['ssid'] for j_ss_info in ss_info]
        final_ss_dict = {'final_ss': final_ss_id}

        # get data set details
        # get estimate id
        estimate_id = [j_ss_info['estimate_id'] for j_ss_info in ss_info]

        # get sample id
        sample_id = [j_ss_info['sample_id'] for j_ss_info in ss_info]

        # get data set id
        data_set_id = [j_ss_info['data_set_id'] for j_ss_info in ss_info]

        # get perturbation id
        perturbation_id = [j_ss_info['perturbation_id'] for j_ss_info in ss_info]

        # collate into single dict
        details_dict = dict(zip(['estimate_id', 'sample_name', 'data_set_id', 'perturbation_id'],
                                [estimate_id, sample_id, data_set_id, perturbation_id]))

        # get everything into single dictionary
        ss_info_dict = {}
        ss_info_dict.update(y_dict)
        ss_info_dict.update(f_dict)
        ss_info_dict.update(final_ss_dict)
        ss_info_dict.update(details_dict)

        self.perturbation_ss = ss_info_dict

        return None

    def __add_initial_ss_id(self, initial_ss):
        """add ss_id from initial wt sims to data dict"""

        wt_ss_id = [j_ss_info['ssid'] for j_data in self.perturbation_ss['estimate_id']
                    for j_ss_info in initial_ss if j_data == j_ss_info['estimate_id']]

        # get ss_id information for initial_ss
        # wt_ss_id = [j_ss_info['ssid'] for j_ss_info in initial_ss]
        wt_ss_dict = {'initial_ss': wt_ss_id}

        self.perturbation_ss.update(wt_ss_dict)

        return None

    def create_ss_perturbation_dict(self, all_results):
        """create complete dictionary of details to be converted to df for writing to file"""

        # separate initial simulation values from perturbation simulation values
        initial_sims, perturbation_sims = self.separate_initial_perturbation(all_results)

        # get initial ss values only
        initial_ss, _ = self.separate_ss_dyn(initial_sims)

        # get perturbation ss values only
        perturbation_ss, _ = self.separate_ss_dyn(perturbation_sims)

        self.__convert_to_ss_dict_for_df(perturbation_ss)

        # add ss_id information for initial_ss
        self.__add_initial_ss_id(initial_ss)

        return None

    def create_df(self, all_results, dummy_arg=[]):
        """create df from validation data information in all results"""
        self.create_ss_perturbation_dict(all_results)

        # convert info to multi index data frame for storage and retrieval
        ss_info = self.perturbation_ss
        df_index_tuple = [(i_estimate, i_sample, i_data_set, i_perturbation)
                          for i_estimate, i_sample, i_data_set, i_perturbation in
                          zip(ss_info['estimate_id'], ss_info['sample_name'], ss_info['data_set_id'],
                              ss_info['perturbation_id'])]
        # df_index_tuples = [(i_value_sample, i_value_exp) for i_value_sample, i_value_exp in
        #                    zip(self.perturbation_ss["sample_name"], self.perturbation_ss["experiment_id"])]
        multi_index_labels = ['estimate_id', 'sample_name', 'data_set_id', 'experiment_id']
        index = pd.MultiIndex.from_tuples(df_index_tuple, names=multi_index_labels)

        # remove unused column headings from dictionary
        del self.perturbation_ss["estimate_id"]
        del self.perturbation_ss["sample_name"]
        del self.perturbation_ss["data_set_id"]
        del self.perturbation_ss["perturbation_id"]

        # create data frame
        all_ss_df = pd.DataFrame(self.perturbation_ss, index=index, columns=self.perturbation_ss.keys())

        # level depth sorting correction
        all_ss_df.sort_index(level=['estimate_id', 'sample_name', 'data_set_id', 'experiment_id'], inplace=True)\

        self.validate_index_label = multi_index_labels

        # write results to file
        if self.validate_file:
            all_ss_df.to_csv(self.validate_file, index_label=multi_index_labels)
            print('\n  Validation Data written to file \n')

        return all_ss_df, multi_index_labels

    def retrieve_validate_df_from_file(self):
        """retrieve validate df from saved file"""
        # read dataframe from csv file
        validate_df = pd.read_csv(self.validate_file, index_col=self.validate_index_label)

        # level depth sorting correction
        validate_df.sort_index(level=['estimate_id', 'sample_name', 'data_set_id', 'experiment_id'], inplace=True)

        # lexicographic sorting of indices in all_ss_df
        # all_ss_df.sort_index(level='estimate_id', inplace=True)
        # all_ss_df.sort_index(level='sample_name', inplace=True)
        # all_ss_df.sort_index(level='data_set_id', inplace=True)
        # all_ss_df.sort_index(level='experiment_id', inplace=True)

        return validate_df

    def retrieve_exp_df_from_file(self):
        """retrieve experimental data from csv file"""
        # read dataframe from csv file
        experiment_df = pd.read_csv(self.original_exp_file, index_col=self.original_index_label)

        # lexographic depth-based ordering of exp df indices
        experiment_df.sort_index(level=['sample_name', 'experiment_id'], inplace=True)
        # experiment_df.sort_index(level='sample_name', inplace=True)
        # experiment_df.sort_index(level='experiment_id', inplace=True)

        return experiment_df

    def process_validation_data(self):
        """get concentration and flux and compare them between original experimental data
        and that obtained from estimated parameter simulations"""

        # retrieve df from file and lex sort by index
        validate_df = self.retrieve_validate_df_from_file()

        # retrieve original set of experiments from file
        experiment_df = self.retrieve_exp_df_from_file()

        # gather concentrations
        y_names, y_values, y_exp_values = self.gather_validation_data(validate_df, experiment_df, 'metabolite')

        # gather fluxes
        f_names, f_values, f_exp_values = self.gather_validation_data(validate_df, experiment_df, 'flux')

        concentration_data = {'names': y_names, 'values': y_values, 'experiment_values': y_exp_values}
        flux_data = {'names': f_names, 'values': f_values, 'experiment_values': f_exp_values}

        self.concentration_validation = concentration_data
        self.flux_validation = flux_data

        return concentration_data, flux_data

    def process_experiment_based_data(self):
        """get concentrations/fluxes for each experiment in each estimate and compare their distributions"""

        # retrieve df from file and lex sort by index
        validate_df = self.retrieve_validate_df_from_file()

        # gather experiment-based validation data
        y_o_names, exp_names, y_o_values = self.ordered_data_collection(validate_df, 'metabolite')

        # gather experiment-based validation data on concentrations
        f_o_names, f_exp_names, f_o_values = self.ordered_data_collection(validate_df, 'flux')

        concentration_data = {'names': y_o_names, 'values': y_o_values, 'experiment_id': exp_names}
        flux_data = {'names': f_o_names, 'values': f_o_values, 'experiment_id': f_exp_names}

        return concentration_data, flux_data

    @staticmethod
    def ordered_data_collection(df, variable_type, select_values=[]):
        """collect concentration/fluxes for each estimate for each experiment"""

        # get variable name
        var_names = variable_name(variable_type)

        # validate_df levels: level[0] - level[3] - estimate, sample, dataset, experiment
        experiment_names = list(df.index.levels[3].unique())
        # number_experiments = len(experiment_names)

        idx = pd.IndexSlice
        df_values = []
        for i_variable in var_names:
            if select_values and i_variable in select_values:
                j_variable_values = []
                for i_experiment in experiment_names:
                    j_variable_values.append(
                        [i_value for i_value in df.loc[idx[:, :, :, i_experiment], i_variable].values])
                df_values.append(j_variable_values)
            else:
                j_variable_values = []
                for i_experiment in experiment_names:
                    j_variable_values.append([i_value
                                              for i_value in df.loc[idx[:, :, :, i_experiment], i_variable].values])
                df_values.append(j_variable_values)

        return var_names, experiment_names, df_values

    @staticmethod
    def gather_validation_data(df, exp_df, variable_type, select_value=[]):
        """gather validation data in orderly fashion for further processing"""

        # get variable name
        var_names = variable_name(variable_type)

        # validate_df levels: level[0] - level[3] - estimate, sample, dataset, experiment
        sample_names = list(df.index.levels[1].unique())
        all_data_sets = df.index.levels[2].unique()
        # number_samples = len(sample_names)
        number_data_sets = len(all_data_sets)

        # experiment_df levels: level[0] - level[1] - sample, experiment
        # exp_sample_names = exp_df.index.levels[0].unique()
        # number_exp_samples = len(exp_sample_names)

        idx = pd.IndexSlice
        if select_value:
            df_values = [[i_value for i_sample in sample_names
                          for i_value in df.loc[idx[:, i_sample, :, :], i_variable].values] for i_variable in var_names
                         if i_variable in select_value]
        else:
            df_values = [[i_value for i_sample in sample_names
                          for i_value in df.loc[idx[:, i_sample, :, :], i_variable].values] for i_variable in var_names]

        desired_exp = []
        for i_variable in var_names:
            if select_value and i_variable in select_value:
                desired_exp.append([i_value for i_sample in sample_names
                                    for i_value in exp_df.loc[idx[i_sample, :], i_variable].values] * number_data_sets)
            else:
                desired_exp.append([i_value for i_sample in sample_names
                                    for i_value in exp_df.loc[idx[i_sample, :], i_variable].values] * number_data_sets)

        return var_names, df_values, desired_exp

    @staticmethod
    def validation_plots(info, scatter=True, box=False, violin=True):
        """validation plots for concentrations/fluxes in info"""

        number_variables = len(info["names"])
        if scatter:
            f1 = plt.figure(figsize=(10, 8), dpi=100, tight_layout=True)
            plot_grid = gridspec.GridSpec(3, number_variables)
        else:
            f1 = plt.figure(figsize=(10, 8), dpi=100, tight_layout=True)
            plot_grid = gridspec.GridSpec(2, number_variables)

        # box plot
        if box:
            box_axis = f1.add_subplot(plot_grid[0, :])
            plot_on_axis_object_box(box_axis, info["values"])
            box_axis.set_xticklabels(info["names"])

            # plot histogram
            validation_hist(info["values"], info["names"], figure_object=f1, grid_objects=plot_grid)

            # scatter plot
            if scatter:
                validation_scatter(info, figure_object=f1, grid_objects=plot_grid)

        # violin plot
        if violin:
            violin_axis = f1.add_subplot(plot_grid[0, :])
            plot_on_axis_object_violin(violin_axis, info["values"])
            violin_axis.set_xticks(np.arange(1, len(info["names"]) + 1))
            violin_axis.set_xticklabels(info["names"])

            # plot histogram
            validation_hist(info["values"], info["names"], figure_object=f1, grid_objects=plot_grid)

            # scatter plot
            if scatter:
                validation_scatter(info, figure_object=f1, grid_objects=plot_grid)

        return None

    def separate_validation_plot(self, scatter=True, box=False, violin=True):
        """plot scatter, hist and box/violin plot for concentration/fluxes determined
        through validation simulations"""

        self.process_validation_data()

        # plot validated concentrations
        self.validation_plots(self.concentration_validation, scatter=scatter, box=box, violin=violin)

        # plot validated fluxes
        self.validation_plots(self.flux_validation, scatter=scatter, box=box, violin=violin)


        # number_variables = len(info_dict["names"])
        # if scatter:
        #     f1 = plt.figure(figsize=(10, 8), dpi=100, tight_layout=True)
        #     plot_grid = gridspec.GridSpec(3, number_variables)
        # else:
        #     f1 = plt.figure(figsize=(10, 8), dpi=100, tight_layout=True)
        #     plot_grid = gridspec.GridSpec(2, number_variables)
        #
        # # box plot
        # if box:
        #     box_axis = f1.add_subplot(plot_grid[0, :])
        #     plot_on_axis_object_box(box_axis, info_dict["values"])
        #     box_axis.set_xticklabels(info_dict["names"])
        #
        #     # plot histogram
        #     validation_hist(info_dict["values"], info_dict["names"], figure_object=f1, grid_objects=plot_grid)
        #
        #     # scatter plot
        #     if scatter:
        #         validation_scatter(info_dict, figure_object=f1, grid_objects=plot_grid)
        #
        # # violin plot
        # if violin:
        #     violin_axis = f1.add_subplot(plot_grid[0, :])
        #     plot_on_axis_object_violin(violin_axis, info_dict["values"])
        #     violin_axis.set_xticks(np.arange(1, len(info_dict["names"]) + 1))
        #     violin_axis.set_xticklabels(info_dict["names"])
        #
        #     # plot histogram
        #     validation_hist(info_dict["values"], info_dict["names"], figure_object=f1, grid_objects=plot_grid)
        #
        #     # scatter plot
        #     if scatter:
        #         validation_scatter(info_dict, figure_object=f1, grid_objects=plot_grid)
        return None
