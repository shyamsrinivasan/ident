import numpy as np
import itertools as it
import pandas as pd
import matplotlib.pyplot as plt
from names_strings import true_parameter_values
from names_strings import ident_parameter_name
from plot_ident_results import plot_on_axis_object_polar


def write_ident_info_file(all_data_dict, exp_df, file_name):
    """create data frame from identifiability data and write to csv file for future use"""
    # reset index of experimental data df
    reset_exp_df = exp_df.reset_index("sample_name")
    reset_exp_df.sort_index(level='data_set_id', inplace=True)
    reset_exp_df.reset_index("data_set_id", inplace=True)
    # reset_exp_df = exp_df.reset_index('experiment_id')
    # reset_exp_df.reset_index('sample_name', inplace=True)
    # lexographic ordering of df indices    #
    # reset_exp_df.sort_index(level='sample_name', inplace=True)

    # create data frame
    data_df = pd.DataFrame(all_data_dict, columns=all_data_dict.keys())

    # number of occurrences of each data set id = number of experiments per data set in the first sample
    first_sample_rows = data_df[data_df["sample_name"] == 'sample_0']
    data_set_id_frequency = int(max(first_sample_rows["data_set_id"].value_counts()))
    # all experiment ids
    experiment_pos_names = ['experiment_{}_id'.format(i_experiment) for i_experiment in range(0, data_set_id_frequency)]
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
                              for i_p in [j_data_set[j_position_exp]]*len(experiment_pos_names)]*5
                             for j_position_exp in range(0, len(experiment_pos_names))]
    all_pos_exp_parameters = [[i_p for j_data_set in all_data_set_exp_parameters
                               for i_p in [j_data_set[j_position_exp]] * len(experiment_pos_parameters)] * 5
                              for j_position_exp in range(0, len(experiment_pos_parameters))]
    experiment_pos_info_keys = experiment_pos_names + experiment_pos_parameters
    experiment_pos_info_values = all_pos_experiment_id + all_pos_exp_parameters
    exp_info_dict = dict(zip(experiment_pos_info_keys, experiment_pos_info_values))
    all_data_dict.update(exp_info_dict)

    # multi index tuples
    ind_tuple = [(j_sample, j_data_set) for j_sample, j_data_set in
                 zip(all_data_dict["sample_name"], all_data_dict["data_set_id"])]

    # multi index index
    index_label = ['sample_name', 'data_set_id']
    index = pd.MultiIndex.from_tuples(ind_tuple, names=index_label)

    # remove redundant columns
    del all_data_dict["sample_name"]
    del all_data_dict["data_set_id"]

    # create multi index data frame
    all_data_df = pd.DataFrame(all_data_dict, index=index, columns=all_data_dict.keys())

    # save data frame to csv file
    all_data_df.to_csv(file_name, index_label=index_label)
    return all_data_df, index_label


def get_all_indices(mother_list, value):
    current_value_id = []
    for i in it.compress(range(0, len(mother_list)),
                         [True if j_id == value else False for j_id in mother_list]):
        current_value_id.append(i)
    return current_value_id


def get_most_useful_dataset(ident_boolean_array):
    """get data identifying most parameters (most useful data set)"""
    # most useful dataset - based on number of parameter identified
    number_data, number_parameter = ident_boolean_array.shape
    identified_parameters = []  # np.zeros((number_data, 1))
    for idata in range(0, number_data):
        ip = [p_id for p_id in it.compress(range(0, number_parameter), ident_boolean_array[idata, :])]
        identified_parameters.append(len(ip))

    # get data identifying most parameters (most useful data set)
    max_useful_data = max(identified_parameters)
    max_data_id = get_all_indices(identified_parameters, max_useful_data)
    max_data = {"maximum": max_useful_data,
                "id": max_data_id,
                "parameters": identified_parameters}
    return max_data


def get_top_useful_data(max_data):
    number_parameters_ided = max_data["parameters"]
    maximum_val = max(number_parameters_ided)
    all_current_values = []
    all_data_id = []
    current_value = maximum_val
    while current_value >= min(number_parameters_ided):
        current_value_id = get_all_indices(number_parameters_ided, current_value)
        all_data_id.append(current_value_id)
        all_current_values.append(current_value)
        current_value -= 1
    return all_current_values, all_data_id


def data_usefulness_percentage(ident_details):
    """get percentage of data identifying n, n-1, n-2, ... number of parameters. n is the maximum number of
    parameters identified by a data set"""
    number_data, _ = ident_details["boolean"].shape
    # most useful dataset - based on number of parameter identified
    max_data = get_most_useful_dataset(ident_details["boolean"])

    # get top n, n-1, n-2, .... parameters identifying data sets
    number_parameters_identified, data_id = get_top_useful_data(max_data)

    # percentage of datasets identifying x number of parameters
    total_number_of_data = []
    percentage_of_total = []
    for i, j_parameter_identified in enumerate(number_parameters_identified):
        number_data_identifying_j = len(data_id[i])
        total_number_of_data.append(number_data_identifying_j)
        percentage_of_total.append(float(number_data_identifying_j) / float(number_data) * 100)
    data_usefulness = {'data set size': number_data,
                       'number': number_parameters_identified,  # number of parameters ided
                       'index': data_id,                        # index of data combinations that id x parameters
                       'total': total_number_of_data,           # number of data combinations that id x parameters
                       'percentage': percentage_of_total,       # percentage of data combinations that id x parameters
                       'number of data combinations': number_data,  # total number of data combinations used
                       'flux id': ident_details["flux id"],         # flux id of data being processed
                       'flux choice': ident_details["flux choice"]} # choice of flux version being identified
    return data_usefulness


def data_utility(ident_details):
    """get information on how useful each data combination is based on the number of parameters in each flux that each
    data combination can identify"""

    # get data identification percentages to classify utility of data sets
    data_usefulness = data_usefulness_percentage(ident_details)

    # print total individual data sets required and combinations found
    number_original_data = 0
    for list_pos, i_data in enumerate(data_usefulness["index"]):
        # original data sets
        print('Original Data sets that can detect {} parameters: {}'.
            format(data_usefulness["number"][list_pos], len(i_data)))
        number_original_data += len(i_data)
    return data_usefulness


def collate_flux_based_data(sample_ident_detail):
    """get processed data for inidividual flux within each individual sample as input args
    and collate them for all parameters in all fluxes within the sample"""
    number_of_fluxes = len(sample_ident_detail)
    for i_flux, i_flux_detail in enumerate(sample_ident_detail):
        # boolean stack
        try:
            all_flux_boolean = np.hstack((all_flux_boolean, i_flux_detail["boolean"]))
        except NameError:
            all_flux_boolean = i_flux_detail["boolean"]
        # nr, dr and final value (identifiability) stack
        try:
            all_flux_values = np.hstack((all_flux_values, i_flux_detail["values"]))
        except NameError:
            all_flux_values = i_flux_detail["values"]
        # number of parameters stack
        try:
            all_flux_parameters = np.hstack((all_flux_parameters, i_flux_detail["parameters"]))
        except NameError:
            all_flux_parameters = i_flux_detail["parameters"]
        # flux id stack
        try:
            flux_id_list.append(i_flux_detail["flux id"])
        except NameError:
            flux_id_list = [i_flux_detail["flux id"]]
        try:
            flux_choice_list.append(i_flux_detail["flux choice"])
        except NameError:
            flux_choice_list = [i_flux_detail["flux choice"]]
    combined_flux_ident_data = {"boolean": all_flux_boolean,
                                "values": all_flux_values,
                                "parameters": all_flux_parameters,
                                "flux id": flux_id_list,
                                "flux choice": flux_choice_list}
    return combined_flux_ident_data


def combined_sample_based_averages_data_utility(all_sample_combined_flux_data_list):
    """sample based averages of data utility when parameters from
    all fluxes are combined into one big pot"""
    j_sample_total = []
    j_sample_percent = []
    j_sample_number = []
    for j_sample_id, j_sample_data in enumerate(all_sample_combined_flux_data_list):
        j_sample_total.append(j_sample_data["total"])
        j_sample_percent.append(j_sample_data["percentage"])
        j_sample_number.append(j_sample_data["number"])
    sample_total_mean = list(np.mean(np.array(j_sample_total), axis=0))
    sample_total_std = list(np.std(np.array(j_sample_total), axis=0))
    sample_percent_mean = list(np.mean(np.array(j_sample_percent), axis=0))
    sample_percent_std = list(np.std(np.array(j_sample_percent), axis=0))
    sample_number_mean = list(np.mean(np.array(j_sample_number), axis=0))
    all_processed_total = {"mean": sample_total_mean,
                           "std": sample_total_std,
                           "number": sample_number_mean,
                           "total data": all_sample_combined_flux_data_list[0]["number of data combinations"],
                           "flux id": all_sample_combined_flux_data_list[0]["flux id"],
                           "flux choice": all_sample_combined_flux_data_list[0]["flux choice"]}
    all_processed_percent = {"mean": sample_percent_mean,
                             "std": sample_percent_std,
                             "number": sample_number_mean,
                             "total data": all_sample_combined_flux_data_list[0]["number of data combinations"],
                             "flux id": all_sample_combined_flux_data_list[0]["flux id"],
                             "flux choice": all_sample_combined_flux_data_list[0]["flux choice"]}
    return {"total": all_processed_total, "percent": all_processed_percent}


def parameter_ident_info(ident_df):
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

    all_p_info = {"names": all_parameter_names,
                  "values": all_p_values,
                  "identifiability": all_p_identifiability,
                  "data_sets": all_p_data_set_names,
                  "total_data_sets": total_data_sets,
                  "identifiability_percentage": [np.array(float(i_nr) * 100/float(i_dr)) for i_nr, i_dr in
                                                 zip(all_p_identifiability, total_data_sets)],
                  "sample_name": all_p_sample_names}
    return all_p_info


def sample_ident_info(all_sample_info):
    """return calculated mean identifiability and std in identifiability
    from all samples for all parameters"""
    # get common data sets identifying all parameters in all samples
    number_parameters = len(all_sample_info[0]["names"])
    # all_sample_data_sets = [[i_sample_info["data_sets"][i_parameter] for i_sample_info in all_sample_info]
    #                         for i_parameter in range(0, number_parameters)]
    # all_parameter_data_sets = []
    # for i_parameter in all_sample_data_sets:
    #     i_parameter_common_data_set = set(i_parameter[0])
    #     for i_sample_info in i_parameter[1:]:
    #         i_parameter_common_data_set.intersection_update(i_sample_info)
    #     all_parameter_data_sets.append(list(i_parameter_common_data_set))

    all_identifiabilites = [i_sample_info["identifiability"] for i_sample_info in all_sample_info]
    sample_mean_identifiability = np.mean(np.array(all_identifiabilites), axis=0)
    sample_mean_ident = [np.array(i_parameter_ident) for i_parameter_ident in sample_mean_identifiability]
    sample_std_identifiability = np.std(np.array(all_identifiabilites), axis=0)
    sample_std_ident = [np.array(i_parameter_ident) for i_parameter_ident in sample_std_identifiability]
    all_identifiabilites_percent = [i_sample_info["identifiability_percentage"]
                                    for i_sample_info in all_sample_info]
    sample_mean_identifiability_percent = np.mean(np.array(all_identifiabilites_percent), axis=0)
    sample_mean_ident_percent = [np.array(i_parameter_ident) for i_parameter_ident in
                                 sample_mean_identifiability_percent]
    sample_std_identifiability_percent = np.std(np.array(all_identifiabilites_percent), axis=0)
    sample_std_ident_percent = [np.array(i_parameter_ident) for i_parameter_ident in
                                sample_std_identifiability_percent]
    all_sample_data_pair = [[i_p for i_sample in all_sample_info
                             for i_p in zip(i_sample["sample_name"][i_parameter], i_sample["data_sets"][i_parameter])]
                            for i_parameter in range(0, number_parameters)]
    ident_dict = {"ident_mean": sample_mean_ident,
                  "ident_std": sample_std_ident,
                  "ident_percent_mean": sample_mean_ident_percent,
                  "ident_percent_std": sample_std_ident_percent,
                  "sample_data_set_id": all_sample_data_pair}
    return ident_dict


def parameter_exp_info(ident_df, exp_df, parameter_ident_info):
    """return experiment type information for each parameter"""
    # get parameter names
    all_parameter_names = ident_df["parameter_name"].unique().tolist()
    number_experiments = len(all_parameter_names)
    exp_column_ids = ['experiment_{}_parameter'.format(i_experiment) for i_experiment in range(0, number_experiments)]
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
            name_value_pair = name_value_pair + zip(missing_perturbation, [0.0] * len(missing_perturbation))
            # names, values = map(list, zip(*name_value_pair))
            # arrange name/values in desired order for every parameter
            given_value = []
            for i_given_name in all_possible_perturbations:
                given_value.append([i_obtained_value for i_obtained_name, i_obtained_value in name_value_pair
                                    if i_given_name == i_obtained_name][0])
            all_experiment_info.update({exp_column_name[i_experiment]: {"names": list(all_possible_perturbations),
                                                                        "frequency": given_value}})
        all_parameter_exp_info.append(all_experiment_info)

    return all_parameter_exp_info


def process_ident(ident_df, exp_df, ident=1, exp_info=1):
    """process ident data to create final data frame of results for plotting"""

    idx = pd.IndexSlice

    # lexographic ordering of exp df indices
    exp_df.sort_index(level='sample_name', inplace=True)
    exp_df.sort_index(level='data_set_id', inplace=True)
    exp_df.sort_index(level='experiment_id', inplace=True)
    # reset_exp_df = exp_df.reset_index('experiment_id')

    # lexicographic ordering of ident df indices
    ident_df.sort_index(level='sample_name', inplace=True)
    ident_df.sort_index(level='data_set_id', inplace=True)
    # reset_ident_df = ident_df.reset_index('data_set_id')

    # number of samples
    sample_names = ident_df.index.levels[0].values.tolist()
    number_samples = len(sample_names)

    if ident:
        # parameter identifiability information for each sample
        all_sample_info = []
        for i_sample in sample_names:
            j_sample_info = ident_df.loc[idx[i_sample, :], :]
            # get all data sets identifying each parameter in each sample
            all_p_info = parameter_ident_info(j_sample_info)
            all_sample_info.append(all_p_info)

        ident_dict = sample_ident_info(all_sample_info)
    else:
        ident_dict = {}
    # parameter identifiability information for all samples along with mean and std between samples
    all_parameter_info = parameter_ident_info(ident_df)
    all_parameter_info.update(ident_dict)

    # get experiment information for each parameter (only for noise-less data)
    # number_experiments = len(all_parameter_info["names"])
    if exp_info and number_samples == 1:
        exp_info = parameter_exp_info(ident_df, exp_df, all_parameter_info["ident_mean"])
    else:
        exp_info = []
    all_parameter_info.update({"exp_info": exp_info})

    # get flux names
    all_flux_names = ident_df["flux_name"].unique().tolist() * len(all_parameter_info["names"])
    all_parameter_info.update({"flux_name": all_flux_names})

    return all_parameter_info


def get_parameter_value(info_dict, ident_df):
    """extract parameter values in a given flux to re-simulate model with newly determined parameters.
    get parameter values from data sets that can detect all parameters"""

    # lexicographic ordering of ident df indices
    ident_df.sort_index(level='sample_name', inplace=True)
    ident_df.sort_index(level='data_set_id', inplace=True)

    # get data sets identifying each parameter
    identifying_data_sets = [set(i_parameter_data_set) for i_parameter_data_set in info_dict["sample_data_set_id"]]
    size_of_data_sets = [len(i_parameter_set) for i_parameter_set in identifying_data_sets]
    sort_index = np.argsort(size_of_data_sets)  # last position is the biggest data set
    largest_set = identifying_data_sets[sort_index[-1]]
    # del identifying_data_sets[sort_index[-1]]
    for i_index in range(len(sort_index)-2, -1, -1):
        largest_set.intersection_update(identifying_data_sets[sort_index[i_index]])

    # get parameter values from data sets in largest_set
    idx = pd.IndexSlice
    relevant_df = ident_df.loc[idx[list(largest_set)], ['parameter_name', 'parameter_value']]
    all_parameter_values = []
    all_parameter_names = []
    all_parameter_data_sets = [list(largest_set) for _ in info_dict["names"]]
    for i_parameter, i_parameter_name in enumerate(info_dict["names"]):
        i_p_value = relevant_df[relevant_df["parameter_name"] == i_parameter_name]["parameter_value"].values
        all_parameter_values.append(i_p_value)
        all_parameter_names.append(i_parameter_name)
    all_parameter_info = {"names": all_parameter_names, "values": all_parameter_values,
                          "data_sets": all_parameter_data_sets, "total_data_sets": info_dict["total_data_sets"],
                          "flux_name": info_dict["flux_name"]}

    return all_parameter_info


def extract_parameter_values_numerical(all_parameter_info):
    """extract parameter values from numerical identifiability method"""
    number_parameters = len(all_parameter_info["values"])
    if number_parameters > 0:
        number_estimates = len(all_parameter_info["values"][0])
    else:
        number_estimates = 0
    parameter_names = all_parameter_info["names"]
    parameter_values = all_parameter_info["values"]
    all_parameter_estimates = [[np.array(i_parameter[j_parameter_estimate]) for i_parameter in parameter_values]
                               for j_parameter_estimate in range(0, number_estimates)]
    extracted_parameter_values = [dict(zip(parameter_names, i_parameter_estimate))
                                  for i_parameter_estimate in all_parameter_estimates]
    extracted_parameters = [{"parameter_name": parameter_names,
                             "parameter_value": extracted_parameter_values,
                             "data_id": all_parameter_info["data_id"]}]
    return extracted_parameters
