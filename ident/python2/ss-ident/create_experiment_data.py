import numpy as np
import pandas as pd
from generate_expdata import generate_expdata
from simulate_data import arrange_experimental_data


def create_experiment_data(save_file_name, noise=0, kinetics=2, number_samples=1, noise_std=0.05):
    """run generate_expdata and store resulting experimental data as csv file from data frame"""
    # generate no-noise experimental data for testing identifiability
    y0 = np.array([5, 1, 1])
    # default parameter values
    cvode_options = ('Newton', 'Adams', 1e-10, 1e-10, 200)
    ode_parameter_values = {"K1ac": np.array([.1]),
                            "K3fdp": np.array([.1]),
                            "L3fdp": np.array([4e6]),
                            "K3pep": np.array([.1]),
                            "K2pep": np.array([.3]),
                            "vemax": np.array([1.1]),
                            "Kefdp": np.array([.45]),
                            "ne": np.array([2]),
                            "d": np.array([.25]),
                            "V4max": np.array([.2]),
                            "k1cat": np.array([1]),
                            "V3max": np.array([1]),
                            "V2max": np.array([1]),
                            "ac": np.array([.1])}

    # get experimental system steady state data without noise using Convenience Kinetics for v3 (kinetics = 2)
    experiment_df, multi_index_labels = generate_expdata(y0, cvode_options, ode_parameter_values,
                                                         noise=noise, kinetics=kinetics,
                                                         dynamic_plot=0, perturbation_plot=0,
                                                         number_of_samples=number_samples, noise_std=noise_std)
    # save data frame to csv file
    experiment_df.to_csv(save_file_name, index_label=multi_index_labels)
    return experiment_df


def retrieve_experimental_data_from_file(data_file_name, multi_index_label):
    """retrieve experimental data from csv file"""
    # read dataframe from csv file
    experiment_df = pd.read_csv(data_file_name, index_col=multi_index_label)
    return experiment_df


def get_variables(data_frame, variable_label):
    """retrieve all concentration/fluxes from data frame"""
    sample_ids = list(data_frame.index.levels[0])
    # experiment_ids = list(data_frame.index.levels[1])
    # number_samples = len(sample_ids)
    # number_experiments = len(experiment_ids)
    all_variables = []
    # sample loop
    for i_sample in sample_ids:
        df_slice = data_frame.loc[(i_sample, slice(None)), variable_label]
        concentration_list = df_slice.values.tolist()
        concentration_array = [np.array(i_concentration_list) for i_concentration_list in concentration_list]
        all_variables.append(concentration_array)

    if 'initial_ss' in variable_label or 'final_ss' in variable_label:
        all_sample_initial_ss = []
        all_sample_final_ss = []
        for i_sample_info in all_variables:
            if 'initial_ss' in variable_label:
                all_initial_ss = [int(i_experiment_info[0]) for i_experiment_info in i_sample_info]
            else:
                all_initial_ss = []
            if 'final_ss' in variable_label:
                all_final_ss = [int(i_experiment_info[1]) for i_experiment_info in i_sample_info]
            else:
                all_final_ss = []
            all_sample_initial_ss.append(all_initial_ss)
            all_sample_final_ss.append(all_final_ss)
        all_variables = {"initial_ss": all_sample_initial_ss, "final_ss": all_sample_final_ss}

    # row_names = list(data_frame.index.names)
    # column_indices = list(data_frame.columns.values)
    # all_concentrations = []
    # all_row_tuples = []
    # for i_label in concentration_label:
    #     try:
    #         all_concentrations.append(data_frame[i_label].values.tolist())
    #         all_row_tuples.append(list(data_frame["pep"].index.values))
    #     except KeyError:
    #         if i_label not in column_indices:
    #             print('Concentration {} not in Data Frame'.format(i_label))
    # # get only one set of row tuples separated into sample ids and experiment ids
    # sample_ids = [i_sample_id for i_sample_id, _ in all_row_tuples[0]]
    # experiment_ids = [i_experiment_id for _, i_experiment_id in all_row_tuples[0]]

    return all_variables


def get_other_details(data_frame, info_label):
    """retrieve miscellaneous details about each data set in data frame (perturbation_details)"""
    sample_ids = list(data_frame.index.levels[0])
    df_slice = data_frame.loc[(sample_ids[0], slice(None)), info_label]
    all_info = df_slice.values.tolist()
    all_info_list = []
    for j_info, _ in enumerate(info_label):
        all_info_list.append([i_experiment_info[j_info] for i_experiment_info in all_info])

    all_info_dict = dict(zip(info_label, all_info_list))
    # all_info_dict["final_ss"] = [int(j_experiment_info) for j_experiment_info in all_info_dict["final_ss"]]
    # all_info_dict["initial_ss"] = [int(j_experiment_info) for j_experiment_info in all_info_dict["initial_ss"]]
    all_info_dict["parameter_change"] = [np.array(j_experiment_info)
                                         for j_experiment_info in all_info_dict["parameter_change"]]
    all_info_dict["parameter_change_percentage"] = [np.array(j_experiment_info) for j_experiment_info in
                                                    all_info_dict["parameter_change_percentage"]]

    return all_info_dict


def extract_info_from_df(data_frame, info_label, sample_name=[], experiment_id=[]):
    """extract whatever info is required from data frame for given column for given experiment for given sample"""
    sample_ids = list(data_frame.index.levels[0])
    experiment_ids = list(data_frame.index.levels[1])

    if sample_name and experiment_id:
        df_slice = data_frame.loc[(sample_name, experiment_id), info_label]
    elif sample_name:
        df_slice = data_frame.loc[(sample_name, slice(None)), info_label]
    elif experiment_id:
        df_slice = data_frame.loc[(slice(None), experiment_id), info_label]
    else:
        df_slice = data_frame.loc[:, info_label]

    # number_samples = len(sample_ids)
    # number_experiments = len(experiment_ids)
    # all_variables = []
    # # sample loop
    # for i_sample in sample_ids:
    #     df_slice = data_frame.loc[(i_sample, slice(None)), variable_label]
    #     concentration_list = df_slice.values.tolist()
    #     concentration_array = [np.array(i_concentration_list) for i_concentration_list in concentration_list]
    #     all_variables.append(concentration_array)
    return df_slice


def retrieve_experimental_data(file_name, multi_index_lablel):
    """retrieve experimental data from csv file and collect concentrations, fluxes and
    other information from resulting data frame and pass as output arguments"""
    experimental_df = retrieve_experimental_data_from_file(file_name, multi_index_lablel)
    # get concentrations
    all_concentrations = get_variables(experimental_df, ['pep', 'fdp', 'E'])
    # get fluxes
    all_fluxes = get_variables(experimental_df, ['v1', 'v2', 'v3', 'v4', 'v5', 'v6'])
    # collate all experimental variable info
    exp_ss = {"y": all_concentrations, "flux": all_fluxes}
    # get perturbation details : chnged parameters, their values, etc
    details_header = ['parameter_name', 'parameter_value', 'parameter_change', 'parameter_change_percentage']
    parameter_info = get_other_details(experimental_df, details_header)
    ss_id_header = ['initial_ss', 'final_ss']
    ss_id_info = get_variables(experimental_df, ss_id_header)
    exp_ss.update(ss_id_info)
    # exp_ss.update(dict(zip(details_header, parameter_info)))
    return exp_ss, parameter_info


def create_data_for_analysis(experimental_data_df, experiment_choice_index, experiments_per_set, data_file_name):
    """arrange experimental data and create data frame and
    store for future use for use with identifiability analysis"""
    all_experiment_indices = ['experiment_0', 'experiment_1', 'experiment_2', 'experiment_3', 'experiment_4',
                              'experiment_5', 'experiment_6', 'experiment_7', 'experiment_8', 'experiment_9',
                              'experiment_10', 'experiment_11', 'experiment_12', 'experiment_13', 'experiment_14',
                              'experiment_15', 'experiment_16', 'experiment_17', 'experiment_18', 'experiment_19',
                              'experiment_20']

    experiment_index_choice = [all_experiment_indices[j_choice] for j_choice in experiment_choice_index]

    # get combinations of experimental datasets
    complete_data_dict, index_labels, experiment_choice = \
        arrange_experimental_data(experimental_data_df, experiments_per_set=experiments_per_set,
                                  experiment_choice=experiment_index_choice)
    # create data frame from data
    # same data frame to file
    df_tuples = [(j_sample, j_data_set, j_experiment) for j_sample, j_data_set, j_experiment in
                 zip(complete_data_dict["sample_name"], complete_data_dict["data_set_id"],
                     complete_data_dict["experiment_id"])]

    # create multi index
    index = pd.MultiIndex.from_tuples(df_tuples, names=index_labels)
    # remove redundant columns
    del complete_data_dict["sample_name"]
    del complete_data_dict["data_set_id"]
    del complete_data_dict["experiment_id"]

    # create data frame
    data_df = pd.DataFrame(complete_data_dict, index=index, columns=complete_data_dict.keys())

    # save data frame to csv file
    data_df.to_csv(data_file_name, index_label=index_labels)
    return data_df


def extract_and_create_data_for_analysis(original_data_file_name, original_index_labels, new_data_file_name,
                                         experiment_choice_index, experiments_per_set):
    """extract existing information from given file and create new file for future analysis"""
    experimental_df = retrieve_experimental_data_from_file(data_file_name=original_data_file_name,
                                                           multi_index_label=original_index_labels)
    create_data_for_analysis(experimental_df, experiment_choice_index, experiments_per_set, new_data_file_name)

    return None


def create_data_for_flux(flux_id, noise, number_samples=1):
    """create and store data sets for each flux separately"""
    if flux_id == 'v1':
        experiment_choice_id = [0,
                                1, 2, 3, 4, 5,
                                11, 12, 13, 14, 15,
                                16, 17, 18, 19, 20]
        multi_index_labels = ['sample_name', 'experiment_id']
        if noise:
            if number_samples == 5:
                file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels' \
                            '\ident\python2\ss-ident\experiments_noise_5_samples'
                new_data_file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModel' \
                                     '\ident\python2\ss-ident\exp_v1_2_experiments_noise_5_samples'
            elif number_samples == 500:
                file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels' \
                            '\ident\python2\ss-ident\experiments_noise_500_samples'
                new_data_file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModel' \
                                     '\ident\python2\ss-ident\exp_v1_2_experiments_noise_500_samples'
        else:
            file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels' \
                        '\ident\python2\ss-ident\experiments'
            new_data_file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModel' \
                                 '\ident\python2\ss-ident\exp_v1_2_experiments'
        experiments_per_set = 2
    elif flux_id == 'v2':
        experiment_choice_id = [0,
                                1, 2, 3, 4, 5,
                                6, 7, 8, 9, 10,
                                11, 12, 13, 14, 15]
        multi_index_labels = ['sample_name', 'experiment_id']
        if noise:
            if number_samples == 5:
                file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels' \
                            '\ident\python2\ss-ident\experiments_noise_5_samples'
                new_data_file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModel' \
                                     '\ident\python2\ss-ident\exp_v2_2_experiments_noise_5_samples'
            elif number_samples == 500:
                file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels' \
                            '\ident\python2\ss-ident\experiments_noise_500_samples'
                new_data_file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModel' \
                                     '\ident\python2\ss-ident\exp_v2_2_experiments_noise_500_samples'
        else:
            file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels' \
                        '\ident\python2\ss-ident\experiments'
            new_data_file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModel' \
                                 '\ident\python2\ss-ident\exp_v2_2_experiments'
        experiments_per_set = 2
    elif flux_id == 'v3':
        experiment_choice_id = [0,
                                1, 2, 3, 4, 5,
                                6, 7, 8, 9, 10,
                                16, 17, 18, 19, 20]
        multi_index_labels = ['sample_name', 'experiment_id']
        if noise:
            if number_samples == 5:
                file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels' \
                            '\ident\python2\ss-ident\experiments_noise_5_samples'
                new_data_file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels'\
                                     '\ident\python2\ss-ident\exp_v3_3_experiments_noise_5_samples'
            elif number_samples == 500:
                file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels' \
                            '\ident\python2\ss-ident\experiments_noise_500_samples'
                new_data_file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels' \
                                     '\ident\python2\ss-ident\exp_v3_3_experiments_noise_500_samples'
        else:
            file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels' \
                        '\ident\python2\ss-ident\experiments'
            new_data_file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels' \
                                 '\ident\python2\ss-ident\exp_v3_3_experiments'
        experiments_per_set = 3
    elif flux_id == 'v5':
        experiment_choice_id = [0,
                                1, 2, 3, 4, 5,
                                6, 7, 8, 9, 10,
                                11, 12, 13, 14, 15,
                                16, 17, 18, 19, 20]
        multi_index_labels = ['sample_name', 'experiment_id']
        if noise:
            if number_samples == 5:
                file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels' \
                            '\ident\python2\ss-ident\experiments_noise_5_samples'
                new_data_file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModel' \
                                     '\ident\python2\ss-ident\exp_v5_2_experiments_noise_5_samples'
            elif number_samples == 500:
                file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels' \
                            '\ident\python2\ss-ident\experiments_noise_500_samples'
                new_data_file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModel' \
                                     '\ident\python2\ss-ident\exp_v5_2_experiments_noise_500_samples'
        else:
            file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels' \
                        '\ident\python2\ss-ident\experiments'
            new_data_file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModel' \
                                 '\ident\python2\ss-ident\exp_v5_2_experiments'
        experiments_per_set = 2

    extract_and_create_data_for_analysis(file_name, multi_index_labels, new_data_file_name,
                                         experiment_choice_id, experiments_per_set=experiments_per_set)
    return None


if __name__ == "__main__":
    file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\ident\python2\ss-ident\experiments'
    experiment_info_df = create_experiment_data(file_name, noise=0, kinetics=2)
    # experiment_info_df = create_experiment_data(file_name, noise=1, kinetics=2, number_samples=500, noise_std=0.05)
