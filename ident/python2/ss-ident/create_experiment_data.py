import numpy as np
import pandas as pd
from generate_expdata import generate_expdata


def create_experiment_data(file_name, noise=0, kinetics=2, number_samples=1, noise_std=0.05):
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
    experiment_df.to_csv(file_name, index_label=multi_index_labels)
    return experiment_df


def retrieve_experimental_data_from_file(file_name, multi_index_label):
    """retrieve experimental data from csv file"""
    # read dataframe from csv file
    experiment_df = pd.read_csv(file_name, index_col=multi_index_label)
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


if __name__ == "__main__":
    file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\ident\python2\ss-ident\experiments_noise_500_samples'
    # experiment_info_df = create_experiment_data(file_name, noise=0, kinetics=2)
    experiment_info_df = create_experiment_data(file_name, noise=1, kinetics=2, number_samples=500, noise_std=0.05)
