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
    reset_exp_df = exp_df.reset_index('experiment_id')
    reset_exp_df.reset_index('sample_name', inplace=True)
    # lexographic ordering of df indices
    reset_exp_df.sort_index(level='data_set_id', inplace=True)
    # reset_exp_df.sort_index(level='sample_name', inplace=True)

    # create data frame
    data_df = pd.DataFrame(all_data_dict, columns=all_data_dict.keys())

    # number of occurrences of each data set id = number of experiments per data set
    data_set_id_frequency = int(max(data_df["data_set_id"].value_counts()))
    # all experiment ids
    experiment_pos_names = ['experiment_{}_id'.format(i_experiment) for i_experiment in range(0, data_set_id_frequency)]
    experiment_pos_parameters = ['experiment_{}_parameter'.format(i_experiment)
                                 for i_experiment in range(0, data_set_id_frequency)]

    # extract experiment ids for each data set
    # get all data set ids
    data_set_ids = reset_exp_df.index.unique().tolist()
    all_data_set_experiments = [reset_exp_df.xs(j_data_set_id)["experiment_id"].values.tolist()
                                for j_data_set_id in data_set_ids]
    all_data_set_exp_parameters = [reset_exp_df.xs(j_data_set_id)["parameter_name"].values.tolist()
                                   for j_data_set_id in data_set_ids]
    # temp = [j_position_exp for j_data_set_exp in all_data_set_experiments for j_position_exp in j_data_set_exp]
    all_pos_experiment_id = [[i_pos_rep_exp for j_data_set_exp in
                              all_data_set_experiments for i_pos_rep_exp in [j_data_set_exp[j_position_exp]]*3]
                             for j_position_exp in range(0, len(experiment_pos_names))]
    all_pos_exp_parameters = [[i_pos_rep_exp for j_data_set_exp in
                              all_data_set_exp_parameters for i_pos_rep_exp in [j_data_set_exp[j_position_exp]]*3]
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


def parameters_ided_by_combination(ident_boolean, data_id):
    """get ids of all parameters ided by given list of data ids based on
    boolean identifiability matrix passed as input args"""
    _, p = ident_boolean.shape
    all_ided_parameters = []
    for j_data_id in data_id:
        ided_parameter_id = [i for i in it.compress(range(0, p), list(ident_boolean[j_data_id, :]))]
        all_ided_parameters.append(ided_parameter_id)
    return all_ided_parameters


def experiments_in_combination(data_id, experiment_details):
    """get experiment ids that are part of any input data combination"""
    all_experiment_ids = []
    for j_data_id in data_id:
        experiment_ids = experiment_details["experiment_id"][j_data_id]
        all_experiment_ids.append(experiment_ids)
    return all_experiment_ids


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


def experiments_in_ident_data(boolean_ident_data, experiment_data, experiment_type_index, flux_id, flux_choice):
    """get all data combinations identifying a given parameter within a given flux (passed as input)
    and identify experiments used within the identifying data combination"""
    number_of_data, number_of_parameter = boolean_ident_data.shape
    number_of_experiment_types = len(experiment_type_index)
    boolean_ident_data = boolean_ident_data.transpose()
    all_parameter_experiment_type_info = []
    for j_parameter, j_parameter_info in enumerate(boolean_ident_data):
        print("Experiment types for parameter {} of {}".format(j_parameter + 1, number_of_parameter))
        # get data combinations identifying parameter j
        data_identifying_parameter_j = [j_data for j_data, bool_value in enumerate(j_parameter_info) if bool_value]

        # get specific id of experiments in data combinations
        experiments_in_identifying_data = [experiment_data["experiments"][j_data_id]
                                           for j_data_id in data_identifying_parameter_j]
        number_of_experiments_in_combination = len(experiment_data["experiments"][0])

        # get unique experiments, experiment type and frequency in each position of combination
        experiments_in_identifying_data = np.array(experiments_in_identifying_data)
        all_position_experiment_info = []
        for j_position_in_combination in range(0, number_of_experiments_in_combination):
            print("Experiment type frequencies in position {} of {}:".format(j_position_in_combination + 1,
                                                                             number_of_experiments_in_combination))
            if experiments_in_identifying_data.size:
                experiment_type_boolean = [[True if exp_id in i_experiment_type else False
                                            for exp_id in experiments_in_identifying_data[:, j_position_in_combination]]
                                           for i_experiment_type in experiment_type_index]
            else:
                experiment_type_boolean = [list(j_parameter_info)] * len(experiment_type_index)

            # get experiment type frequency
            experiment_type_frequency = [sum(i_experiment_type_boolean)
                                         for i_experiment_type_boolean in experiment_type_boolean]
            total_identifying_data = sum(experiment_type_frequency)
            if total_identifying_data:
                experiment_type_percentage = [float(i_frequency)*100/total_identifying_data
                                              for i_frequency in experiment_type_frequency]
            else:
                experiment_type_percentage = [0] * len(experiment_type_frequency)

            # get corresponding data combination ids
            experiment_type_data_id = [[data_identifying_parameter_j[id]
                                        for id, val in enumerate(i_experiment_type_boolean) if val]
                                       for i_experiment_type_boolean in experiment_type_boolean]

            for i_experiment_type in range(0, number_of_experiment_types):
                print("Frequency for type {} of {}: {} of {} data combinations".
                      format(i_experiment_type, number_of_experiment_types,
                             experiment_type_frequency[i_experiment_type], total_identifying_data))

            all_position_experiment_info.append({"boolean": experiment_type_boolean,
                                                 "frequency": experiment_type_frequency,
                                                 "percentage": experiment_type_percentage,
                                                 "data id": experiment_type_data_id,
                                                 "position": j_position_in_combination,
                                                 "flux id": flux_id,
                                                 "flux choice": flux_choice})
        all_parameter_experiment_type_info.append(all_position_experiment_info)
        print("Experiment type analysis for parameter {} complete \n".format(j_parameter + 1))

    return all_parameter_experiment_type_info


def flux_based_experiment_info(sample_ident_info, experiment_details, experiment_type_indices):
    """parse information by looping through each flux present within each sample (passed as input)
    to get experiments identifying each para,eter within each flux"""
    number_of_fluxes = len(sample_ident_info)
    all_flux_experiment_type_info = []
    for j_flux, j_flux_detail in enumerate(sample_ident_info):
        print("Getting experiments identifying parameters in flux {} of {}".format(j_flux + 1, number_of_fluxes))
        parameter_experiment_type_info = experiments_in_ident_data(j_flux_detail["boolean"],
                                                                   experiment_details,
                                                                   experiment_type_indices,
                                                                   flux_id=j_flux_detail["flux id"],
                                                                   flux_choice=j_flux_detail["flux choice"])
        all_flux_experiment_type_info.append(parameter_experiment_type_info)
        print("Identifying experiments for flux {} of {} complete \n".format(j_flux + 1, number_of_fluxes))
    return all_flux_experiment_type_info


def get_data_info(data_list, ident_details, experiment_details):
    """get information on every single data set obtained from data_usefulness_percentage"""
    original_data = []
    for i_data in range(0, len(data_list["number"])):
        # get parameters identified by given data ids in i_data
        all_ided_parameters = parameters_ided_by_combination(ident_details["boolean"],
                                                             data_list["index"][i_data])
        # get experiment ids used in/that form data combination
        all_experiment_ids = experiments_in_combination(data_list["index"][i_data], experiment_details)
        all_data_info = {"parameters": all_ided_parameters,
                         "experiments": all_experiment_ids}
        original_data.append(all_data_info)
    return original_data


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


def parameter_identifiability(ident_details):
    """get information on how identifiable each parameter within a flux is based on
    number of data combinatons that identify the parameter - degree of identifiability of each parameter"""
    number_data, p = ident_details["boolean"].shape
    identifying_data = []
    identifying_data_percentage = []
    for iparameter in range(0, p):
        idata = [data_id for data_id in it.compress(range(0, number_data), ident_details["boolean"][:, iparameter])]
        identifying_data.append(len(idata))
        identifying_data_percentage.append(float(len(idata)) / number_data * 100)

    # get parameters identified by most data sets (most identifiable parameter)
    max_identified = max(identifying_data)
    max_identifying_data_percentage = max(identifying_data_percentage)
    max_parameter_id = [identifying_data.index(max(identifying_data))]
    for iparameter, number_identifying_data in enumerate(identifying_data):
        if number_identifying_data == max_identified and iparameter not in max_parameter_id:
            max_parameter_id.append(iparameter)
        elif number_identifying_data > max_identified:
            max_identified = number_identifying_data
            max_parameter_id = [iparameter]
    max_parameter = {"maximum number": max_identified,
                     "maximum percentage": max_identifying_data_percentage,
                     "maximum id": max_parameter_id,
                     "info": identifying_data,
                     "percentage": identifying_data_percentage,
                     "flux id": ident_details["flux id"],
                     "flux choice": ident_details["flux choice"],
                     "data set size": number_data}
    return max_parameter


def true_parameter_value(ident_details):
    number_data, _ = ident_details["boolean"].shape
    number_parameters = ident_details["parameters"]
    all_parameter_true_values = []
    # add loop for fluxes  - needed when using combined data from all fluxes
    try:
        total_parameters = 0
        for i_flux in range(0, len(number_parameters)):
            for i_parameter in range(0, number_parameters[i_flux]):
                # get data combinations identifying parameter j
                i_parameter_ident_data = ident_details["boolean"][:, total_parameters]
                data_identifying_parameter_i = [j_data for j_data, bool_value in enumerate(i_parameter_ident_data)
                                                if bool_value]
                # get true values only for identifiable data
                found_value = ident_details["values"][data_identifying_parameter_i, total_parameters, 2]
                # get flux name
                flux_name = 'flux{}'.format(ident_details["flux id"][i_flux])
                flux_choice_id = ident_details["flux choice"][i_flux]
                # get parameter name
                parameter_name = ident_parameter_name(i_parameter, flux_name=flux_name,
                                                      flux_choice_id=flux_choice_id)
                # get true parameter values
                true_value = true_parameter_values(flux_based=1, flux_name=flux_name,
                                                         flux_choice_id=flux_choice_id,
                                                         parameter_id=parameter_name)
                parameter_true_values = {"flux id": ident_details["flux id"],
                                         "flux name": flux_name,
                                         "flux choice": flux_choice_id,
                                         "parameter id": i_parameter,
                                         "parameter name": parameter_name,
                                         "found values": found_value,
                                         "true values": np.tile(true_value, found_value.shape),
                                         "data set size": number_data}
                all_parameter_true_values.append(parameter_true_values)
                total_parameters += 1
    except TypeError:  # when number_parameters is integer and not list
        # eliminate flux-based loop
        for i_parameter in range(0, number_parameters):
            # get data combinations identifying parameter j
            i_parameter_ident_data = ident_details["boolean"][:, i_parameter]
            data_identifying_parameter_i = [j_data for j_data, bool_value in enumerate(i_parameter_ident_data)
                                            if bool_value]
            # get true values only for identifiable data
            found_value = ident_details["values"][data_identifying_parameter_i, i_parameter, 2]
            # get flux name
            try:
                flux_name = ['flux{}'.format(i_flux_id) for i_flux_id in ident_details["flux id"]]
            except TypeError:
                flux_name = 'flux{}'.format(ident_details["flux id"])
            try:
                flux_choice_id = [i_flux_choice_id for i_flux_choice_id in ident_details["flux choice"]]
            except TypeError:
                flux_choice_id = ident_details["flux choice"]
            # get parameter name
            parameter_name = ident_parameter_name(i_parameter, flux_name=flux_name,
                                                  flux_choice_id=flux_choice_id)
            # get true parameter values
            true_value = true_parameter_values(flux_based=1, flux_name=flux_name,
                                                     flux_choice_id=flux_choice_id,
                                                     parameter_id=parameter_name)
            if found_value.size != 0:
                tiled_found_value = np.tile(true_value, found_value.shape)
            else:
                tiled_found_value = true_value
            parameter_true_values = {"flux id": ident_details["flux id"],
                                     "flux name": flux_name,
                                     "flux choice": flux_choice_id,
                                     "parameter id": i_parameter,
                                     "parameter name": parameter_name,
                                     "data id": data_identifying_parameter_i,
                                     "found values": found_value,
                                     "true values": tiled_found_value,
                                     "data set size": number_data}
            all_parameter_true_values.append(parameter_true_values)

    return all_parameter_true_values


def process_info(ident_details):
    """get data utility and parameter identifiability and other data for
    parameters of each flux passed as input"""
    # get data identification percentages to classify utility of data sets
    data_list = data_utility(ident_details)
    # most easily identifiable parameter - based on frequency of identification
    max_parameter = parameter_identifiability(ident_details)
    # get true parameter values (the value of g = Nr/Dr)
    parameter_true_values = true_parameter_value(ident_details)

    return data_list, max_parameter, parameter_true_values


def flux_based_ident_info(sample_ident_detail):
    """parse information by looping though each flux present within each sample that is passed an input argument"""
    number_of_fluxes = len(sample_ident_detail)
    all_flux_data_list = []
    all_flux_max_parameter = []
    all_flux_true_parameter = []
    for j_flux, j_flux_info in enumerate(sample_ident_detail):
        print("Processing identifiability for flux {} of {}".format(j_flux + 1, number_of_fluxes))
        data_list, max_parameter, parameter_true_value = process_info(j_flux_info)
        all_flux_data_list.append(data_list)
        all_flux_max_parameter.append(max_parameter)
        all_flux_true_parameter.append(parameter_true_value)
        print("Information Processing Complete for flux {} \n".format(j_flux + 1))
    return all_flux_data_list, all_flux_max_parameter, all_flux_true_parameter


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


def collate_sample_based_data_utility(number_of_fluxes_per_sample, all_sample_data_list):
    """collate data utility information based on each flux from multiple samples and generate
    averages and standard deviations for each flux"""
    all_flux_data_list = []
    for j_flux in range(0, number_of_fluxes_per_sample[0]):
        j_flux_data_number = []
        j_flux_data_total = []
        j_flux_data_percentage = []
        j_flux_data_indices = []
        j_flux_data_combination_total = []
        for j_sample_id, j_sample_data in enumerate(all_sample_data_list):
            j_flux_info = j_sample_data[j_flux]
            j_flux_data_number.append(j_flux_info["number"])
            j_flux_data_total.append(j_flux_info["total"])
            j_flux_data_percentage.append(j_flux_info["percentage"])
            j_flux_data_combination_total.append(j_flux_info["number of data combinations"])
            j_flux_data_indices.append(j_flux_info["index"])
        j_flux_total_mean = list(np.mean(np.array(j_flux_data_total), axis=0))
        j_flux_total_std = list(np.std(np.array(j_flux_data_total), axis=0))
        j_flux_percent_mean = list(np.mean(np.array(j_flux_data_percentage), axis=0))
        j_flux_percent_std = list(np.std(np.array(j_flux_data_percentage), axis=0))
        j_flux_number_mean = list(np.mean(np.array(j_flux_data_number), axis=0))
        j_flux_processed_total = {"mean": j_flux_total_mean,
                                  "std": j_flux_total_std,
                                  "number": j_flux_number_mean,
                                  "flux id": all_sample_data_list[0][j_flux]["flux id"],
                                  "flux choice": all_sample_data_list[0][j_flux]["flux choice"],
                                  "data set size": all_sample_data_list[0][j_flux]["data set size"]}
        j_flux_processed_percent = {"mean": j_flux_percent_mean,
                                    "std": j_flux_percent_std,
                                    "number": j_flux_number_mean,
                                    "flux id": all_sample_data_list[0][j_flux]["flux id"],
                                    "flux choice": all_sample_data_list[0][j_flux]["flux choice"],
                                    "data set size": all_sample_data_list[0][j_flux]["data set size"]}
        all_flux_data_list.append({"total": j_flux_processed_total,
                                   "percentage": j_flux_processed_percent})
    return all_flux_data_list


def collate_sample_based_identifibaility(number_of_fluxes_per_sample, all_sample_max_parameter):
    """collect parameter identifiability information for each flux from multiple samples and generate
    averages and standard deviations for each parameter for each flux"""
    all_flux_max_parameter = []
    for j_flux in range(0, number_of_fluxes_per_sample[0]):
        j_flux_data_info = []
        j_flux_data_percent = []
        j_flux_data_maximum_number = []
        j_flux_data_maximum_percent = []
        j_flux_data_maximum_id = []
        for j_sample_id, j_sample_data in enumerate(all_sample_max_parameter):
            j_flux_content = j_sample_data[j_flux]
            j_flux_data_info.append(j_flux_content["info"])
            j_flux_data_percent.append(j_flux_content["percentage"])
            j_flux_data_maximum_number.append(j_flux_content["maximum number"])
            j_flux_data_maximum_percent.append(j_flux_content["maximum percentage"])
            j_flux_data_maximum_id.append(j_flux_content["maximum id"])
        j_flux_info_mean = list(np.mean(np.array(j_flux_data_info), axis=0))
        j_flux_info_std = list(np.std(np.array(j_flux_data_info), axis=0))
        j_flux_percent_mean = list(np.mean(np.array(j_flux_data_percent), axis=0))
        j_flux_percent_std = list(np.std(np.array(j_flux_data_percent), axis=0))
        j_flux_processed_total = {"mean": j_flux_info_mean,
                                  "std": j_flux_info_std,
                                  "flux id": all_sample_max_parameter[0][j_flux]["flux id"],
                                  "flux choice": all_sample_max_parameter[0][j_flux]["flux choice"],
                                  "data set size": all_sample_max_parameter[0][j_flux]["data set size"]}
        j_flux_processed_percent = {"mean": j_flux_percent_mean,
                                    "std": j_flux_percent_std,
                                    "flux id": all_sample_max_parameter[0][j_flux]["flux id"],
                                    "flux choice": all_sample_max_parameter[0][j_flux]["flux choice"],
                                    "data set size": all_sample_max_parameter[0][j_flux]["data set size"]}
        all_flux_max_parameter.append({"total": j_flux_processed_total,
                                       "percentage": j_flux_processed_percent})
    return all_flux_max_parameter


def collate_sample_based_experiment_info(number_of_fluxes_per_sample, all_sample_experiment_list):
    """collate information on experiments appearing at each position in
    a data combination for each parameter for each flux from all samples"""
    all_flux_experiment_info = []
    # loop through flux
    for j_flux in range(0, number_of_fluxes_per_sample[0]):
        number_of_parameters_per_flux = len(all_sample_experiment_list[0][j_flux])
        number_of_experiments_used_per_flux = len(all_sample_experiment_list[0][j_flux][0])
        all_parameter_all_position_info = []
        # loop through parameter of each flux
        for k_parameter in range(0, number_of_parameters_per_flux):
            all_position_info = []
            # loop through each experiment position in each parameter in each flux
            for i_position in range(0, number_of_experiments_used_per_flux):
                i_position_frequency = []
                i_position_percent = []
                i_position_data_id = []
                # loop through each sample of noisy data to collate data on each position
                # for each parameter for each flux for all samples
                for j_sample_id, j_sample_data in enumerate(all_sample_experiment_list):
                    i_position_frequency.append(j_sample_data[j_flux][k_parameter][i_position]["frequency"])
                    i_position_percent.append(j_sample_data[j_flux][k_parameter][i_position]["percentage"])
                    i_position_data_id.append(j_sample_data[j_flux][k_parameter][i_position]["data id"])
                # caluclate means and standard deviation of experiment frequency
                # appearing at each position for each parameter for each flux
                i_position_frequency_mean = list(np.mean(np.array(i_position_frequency), axis=0))
                i_position_frequency_std = list(np.std(np.array(i_position_frequency), axis=0))
                i_position_percent_mean = list(np.mean(np.array(i_position_percent), axis=0))
                i_position_percent_std = list(np.std(np.array(i_position_percent), axis=0))
                i_position_frequency_info = {"mean": i_position_frequency_mean,
                                             "std": i_position_frequency_std,
                                             "flux id": all_sample_experiment_list[0][j_flux]
                                             [k_parameter][i_position]["flux id"],
                                             "flux choice": all_sample_experiment_list[0][j_flux]
                                             [k_parameter][i_position]["flux choice"]}
                i_position_percent_info = {"mean": i_position_percent_mean,
                                           "std": i_position_percent_std,
                                           "flux id": all_sample_experiment_list[0][j_flux]
                                           [k_parameter][i_position]["flux id"],
                                           "flux choice": all_sample_experiment_list[0][j_flux]
                                           [k_parameter][i_position]["flux choice"]}
                # collect data on each position
                all_position_info.append({"total": i_position_frequency_info,
                                          "percentage": i_position_percent_info,
                                          "position": i_position})
            # collect info on each parameter in all positions
            all_parameter_all_position_info.append(all_position_info)
        # collect data on all parameters for each flux
        all_flux_experiment_info.append(all_parameter_all_position_info)

    return all_flux_experiment_info


def collate_sample_based_parameter_value(number_of_fluxes_per_sample, all_sample_parameter_value):
    number_samples = len(all_sample_parameter_value)
    if number_samples > 1:
        # if more than one sample calculate mean based on common ids that can
        # identify each parameter between different samples
        all_flux_parameter_value = []
        for j_flux in range(0, number_of_fluxes_per_sample[0]):
            number_of_parameters_per_flux = len(all_sample_parameter_value[0][j_flux])
            all_parameter_info = []
            for k_parameter in range(0, number_of_parameters_per_flux):
                data_set_size = all_sample_parameter_value[0][j_flux][k_parameter]["data set size"]
                all_sample_data_ids = [j_sample_data[j_flux][k_parameter]["data id"]
                                       for j_sample_data in all_sample_parameter_value]
                all_sample_ident_value = [j_sample_ident[j_flux][k_parameter]["found values"]
                                          for j_sample_ident in all_sample_parameter_value]
                k_parameter_true_values = [np.mean(j_sample_ident[j_flux][k_parameter]["true values"])
                                           for j_sample_ident in all_sample_parameter_value]
                k_parameter_true_array = np.array(k_parameter_true_values)
                all_sample_data_id_ident_value = [zip(j_sample_data_ids, j_sample_ident_value)
                                                  for j_sample_data_ids, j_sample_ident_value in
                                                  zip(all_sample_data_ids, all_sample_ident_value)]
                all_sample_boolean = [[True if j_data_id in set(j_sample_data_id) else False
                                       for j_data_id in range(0, data_set_size)]
                                      for j_sample_data_id in all_sample_data_ids]
                data_sample_ident = np.zeros((number_samples, data_set_size))
                for j_sample_id, j_sample_data in enumerate(all_sample_data_id_ident_value):
                    for j_ident_data in j_sample_data:
                        data_id, ident_data = j_ident_data
                        data_sample_ident[j_sample_id, data_id] = ident_data
                # box plot of across sample variations for each data
                all_boolean_sum = sum(np.array(all_sample_boolean))
                bool_accept = [True if value == number_samples else False for value in all_boolean_sum]
                # import matplotlib.pyplot as plt
                # f, ax = plt.subplots(1, 1)
                # ax.boxplot(data_sample_ident[:, bool_accept])
                mean_across_samples = np.mean(data_sample_ident[:, bool_accept], axis=0)
                std_across_samples = np.std(data_sample_ident[:, bool_accept], axis=0)
                mean_across_data = np.mean(data_sample_ident[:, bool_accept], axis=1)
                std_across_data = np.std(data_sample_ident[:, bool_accept], axis=1)
                k_parameter_info = {"sample mean": mean_across_samples,
                                    "sample std": std_across_samples,
                                    "data mean": mean_across_data,
                                    "data std": std_across_data,
                                    "raw sample ident data": data_sample_ident,
                                    "raw sample boolean data": all_sample_boolean,
                                    "accepted sample ident data": data_sample_ident[:, bool_accept],
                                    "data sample mean": [],
                                    "data sample std": [],
                                    "sample data mean": [],
                                    "sample data std": [],
                                    "parameter id": all_sample_parameter_value[0][j_flux][k_parameter]["parameter id"],
                                    "parameter name": all_sample_parameter_value[0][j_flux][k_parameter]["parameter name"],
                                    "flux name": all_sample_parameter_value[0][j_flux][k_parameter]["flux name"],
                                    "flux id": all_sample_parameter_value[0][j_flux][k_parameter]["flux id"],
                                    "flux choice": all_sample_parameter_value[0][j_flux][k_parameter]["flux choice"],
                                    "true value": np.mean(k_parameter_true_array, axis=0)}
                all_parameter_info.append(k_parameter_info)
            all_flux_parameter_value.append(all_parameter_info)
    else:
        all_flux_parameter_value = []
        for j_flux in range(0, number_of_fluxes_per_sample[0]):
            number_of_parameters_per_flux = len(all_sample_parameter_value[0][j_flux])
            all_parameter_info = []
            for k_parameter in range(0, number_of_parameters_per_flux):
                data_set_size = all_sample_parameter_value[0][j_flux][k_parameter]["data set size"]
                all_sample_data_ids = [j_sample_data[j_flux][k_parameter]["data id"]
                                       for j_sample_data in all_sample_parameter_value]
                all_sample_ident_value = [j_sample_ident[j_flux][k_parameter]["found values"]
                                          for j_sample_ident in all_sample_parameter_value]
                all_sample_data_id_ident_value = [zip(j_sample_data_ids, j_sample_ident_value)
                                                  for j_sample_data_ids, j_sample_ident_value in
                                                  zip(all_sample_data_ids, all_sample_ident_value)]
                all_sample_boolean = [[True if j_data_id in set(j_sample_data_id) else False
                                       for j_data_id in range(0, data_set_size)]
                                      for j_sample_data_id in all_sample_data_ids]
                data_sample_ident = np.zeros((number_samples, data_set_size))
                for j_sample_id, j_sample_data in enumerate(all_sample_data_id_ident_value):
                    for j_ident_data in j_sample_data:
                        data_id, ident_data = j_ident_data
                        data_sample_ident[j_sample_id, data_id] = ident_data
                # box plot of across sample variations for each data
                all_boolean_sum = sum(np.array(all_sample_boolean))
                bool_accept = [True if value == number_samples else False for value in all_boolean_sum]
                # get mean across all samples for each data set identifying parameter k
                mean_across_samples = np.mean(data_sample_ident[:, bool_accept], axis=0)
                std_across_samples = np.std(data_sample_ident[:, bool_accept], axis=0)
                mean_across_data = np.mean(data_sample_ident[:, bool_accept], axis=1)
                std_across_data = np.std(data_sample_ident[:, bool_accept], axis=1)
                if data_sample_ident[:, bool_accept].size == 0:
                    mean_across_data = np.array([0.0])
                    std_across_data = np.array([0.0])
                k_parameter_true_values = [j_sample_data[j_flux][k_parameter]["true values"]
                                           for j_sample_data in all_sample_parameter_value]
                # calculate means and standard deviations for each parameter between samples and between data points
                k_parameter_true_array = np.array(k_parameter_true_values)
                # get mean across all data points across all samples
                # mean_across_data_across_samples = np.mean(mean_across_data, axis=0)
                # std_across_data_across_samples = np.std(mean_across_data, axis=0)
                # # get mean across all samples across all data points
                # if mean_across_samples.size != 0:
                #     mean_across_samples_across_data = np.mean(mean_across_samples, axis=0)
                #     std_across_samples_across_data = np.std(mean_across_samples, axis=0)
                # else:
                #     mean_across_samples_across_data = np.array([0.0])
                #     std_across_samples_across_data = np.array([0.0])
                k_parameter_info = {"sample mean": mean_across_samples,
                                    "sample std": std_across_samples,
                                    "data mean": mean_across_data,
                                    "data std": std_across_data,
                                    "raw sample ident data": data_sample_ident,
                                    "raw sample boolean data": all_sample_boolean,
                                    "data sample mean": [],  # mean_across_data_across_samples,
                                    "data sample std": [],  # std_across_data_across_samples,
                                    "sample data mean": [],  # mean_across_samples_across_data,
                                    "sample data std": [],  # std_across_samples_across_data,
                                    "parameter id": all_sample_parameter_value[0][j_flux][k_parameter]["parameter id"],
                                    "parameter name": all_sample_parameter_value[0][j_flux][k_parameter]["parameter name"],
                                    "flux name": all_sample_parameter_value[0][j_flux][k_parameter]["flux name"],
                                    "flux id": all_sample_parameter_value[0][j_flux][k_parameter]["flux id"],
                                    "flux choice": all_sample_parameter_value[0][j_flux][k_parameter]["flux choice"],
                                    "true value": np.mean(np.mean(k_parameter_true_array, axis=1), axis=0)}
                all_parameter_info.append(k_parameter_info)
            all_flux_parameter_value.append(all_parameter_info)
    return all_flux_parameter_value


def sample_based_averages(number_of_fluxes_per_sample,
                          all_sample_data_list, all_sample_max_parameter,
                          all_sample_experiment_list, all_sample_parameter_value):
    """call sample based collate functions and generate averages and standard deviations for
    both data utility and parameter identifiability for each flux"""
    all_flux_data_utility = collate_sample_based_data_utility(number_of_fluxes_per_sample,
                                                              all_sample_data_list)
    all_flux_max_parameter = collate_sample_based_identifibaility(number_of_fluxes_per_sample,
                                                                  all_sample_max_parameter)
    all_flux_experiment_info = collate_sample_based_experiment_info(number_of_fluxes_per_sample,
                                                                    all_sample_experiment_list)
    all_flux_parameter_value = collate_sample_based_parameter_value(number_of_fluxes_per_sample,
                                                                    all_sample_parameter_value)
    return all_flux_data_utility, all_flux_max_parameter, all_flux_experiment_info, all_flux_parameter_value


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


def combined_sample_based_averages_experiment_info(all_sample_combined_flux_experiment_info):
    """process experiment info for multiple samples when
    all fluxes using the same combination of n experimental data"""
    number_of_parameters = len(all_sample_combined_flux_experiment_info[0])
    # loop through parameters
    all_parameter_all_position_info = []
    for j_parameter in range(0, number_of_parameters):
        number_of_experiments_per_parameter = len(all_sample_combined_flux_experiment_info[0][j_parameter])
        # loop through position in data combination (experiment position)
        all_position_info = []
        for i_position in range(0, number_of_experiments_per_parameter):
            i_position_frequency = []
            i_position_percent = []
            for j_sample_id, j_sample_data in enumerate(all_sample_combined_flux_experiment_info):
                i_position_frequency.append(j_sample_data[j_parameter][i_position]["frequency"])
                i_position_percent.append(j_sample_data[j_parameter][i_position]["percentage"])
            i_position_frequency_mean = list(np.mean(np.array(i_position_frequency), axis=0))
            i_position_frequency_std = list(np.std(np.array(i_position_frequency), axis=0))
            i_position_percent_mean = list(np.mean(np.array(i_position_percent), axis=0))
            i_position_percent_std = list(np.std(np.array(i_position_percent), axis=0))
            i_position_frequency_info = {"mean": i_position_frequency_mean,
                                         "std": i_position_frequency_std,
                                         "flux id": all_sample_combined_flux_experiment_info[0]
                                         [j_parameter][i_position]["flux id"],
                                         "flux choice": all_sample_combined_flux_experiment_info[0]
                                         [j_parameter][i_position]["flux choice"]}
            i_position_percent_info = {"mean": i_position_percent_mean,
                                       "std": i_position_percent_std,
                                       "flux id": all_sample_combined_flux_experiment_info[0]
                                       [j_parameter][i_position]["flux id"],
                                       "flux choice": all_sample_combined_flux_experiment_info[0]
                                       [j_parameter][i_position]["flux choice"]}
            # collect data on each position
            all_position_info.append({"frequency": i_position_frequency_info,
                                      "percentage": i_position_percent_info,
                                      "position": i_position})
        # collect data on each parameter
        all_parameter_all_position_info.append(all_position_info)

    return all_parameter_all_position_info


def combined_sampled_based_averages_parameter_value(all_sample_combined_flux_parameter_value):
    """sample based averages for all parameter values when all fluxes using same
    number of data combinations are combined"""
    number_of_parameters = len(all_sample_combined_flux_parameter_value[0])
    all_parameter_info = []
    # loop through parameters
    for j_parameter in range(0, number_of_parameters):
        # loop through samples
        j_parameter_found_values = []
        j_parameter_true_values = []
        for j_sample_id, j_sample_info in enumerate(all_sample_combined_flux_parameter_value):
            j_parameter_found_values.append(j_sample_info[j_parameter]["found values"])
            j_parameter_true_values.append(j_sample_info[j_parameter]["true values"])
        # calculate means and standard deviations for each parameter between samples and between data points
        j_parameter_true_array = np.array(j_parameter_true_values)
        # get mean across all samples for each data set identifying parameter k
        mean_across_samples = np.mean(np.array(j_parameter_found_values), axis=0)
        std_across_samples = np.std(np.array(j_parameter_found_values), axis=0)
        # get mean across al data identifying parameter k for each sample
        mean_across_data = np.mean(np.array(j_parameter_found_values), axis=1)
        std_across_data = np.std(np.array(j_parameter_found_values), axis=1)
        # get mean across all data points across all samples
        mean_across_data_across_samples = np.mean(mean_across_data, axis=0)
        std_across_data_across_samples = np.std(mean_across_data, axis=0)
        # get mean across all samples across all data points
        mean_across_samples_across_data = np.mean(mean_across_samples, axis=0)
        std_across_samples_across_data = np.std(mean_across_samples, axis=0)
        k_parameter_info = {"found values": np.array(j_parameter_found_values),
                            "sample mean": mean_across_samples,
                            "sample std": std_across_samples,
                            "data mean": mean_across_data,
                            "data std": std_across_data,
                            "data sample mean": mean_across_data_across_samples,
                            "data sample std": std_across_data_across_samples,
                            "sample data mean": mean_across_samples_across_data,
                            "sample data std": std_across_samples_across_data,
                            "parameter id": all_sample_combined_flux_parameter_value[0][j_parameter]["parameter id"],
                            "parameter name": all_sample_combined_flux_parameter_value[0][j_parameter]["parameter name"],
                            "flux name": all_sample_combined_flux_parameter_value[0][j_parameter]["flux name"],
                            "flux id": all_sample_combined_flux_parameter_value[0][j_parameter]["flux id"],
                            "flux choice": all_sample_combined_flux_parameter_value[0][j_parameter]["flux choice"],
                            "true value": np.mean(np.mean(j_parameter_true_array, axis=1), axis=0)}
        all_parameter_info.append(k_parameter_info)

    return all_parameter_info


def parameter_ident_info(ident_df):
    """returns info on identifiability and value of each parameter in df"""
    # get parameter names
    all_parameter_names = ident_df["parameter_name"].unique().tolist()

    all_p_values = []
    all_p_data_set_names = []
    all_p_identifiability = []
    total_data_sets = []
    for i_parameter_name in all_parameter_names:
        # get all data sets identifying each parameter
        identifying_df = ident_df[(ident_df["parameter_name"] == i_parameter_name) & (ident_df["identified"])]
        # get data set names
        all_p_data_set_names.append([i_value[1] for i_value in identifying_df.index.values])
        all_p_identifiability.append(len(all_p_data_set_names[-1]))
        # get parameter values
        all_p_values.append([np.array(i_value) for i_value in identifying_df["parameter_value"].values])
        total_data_sets.append(len(ident_df.index.levels[1]))

    all_p_info = {"names": all_parameter_names,
                  "values": all_p_values,
                  "identifiability": all_p_identifiability,
                  "data_sets": all_p_data_set_names,
                  "total_data_sets": total_data_sets,
                  "identifiability_percentage": [np.array(float(i_nr) * 100/float(i_dr)) for i_nr, i_dr in
                                                 zip(all_p_identifiability, total_data_sets)]}
    return all_p_info


def sample_ident_info(all_sample_info):
    """return calculated mean identifiability and std in identifiability
    from all samples for all parameters"""
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
    ident_dict = {"ident_mean": sample_mean_ident,
                  "ident_std": sample_std_ident,
                  "ident_percent_mean": sample_mean_ident_percent,
                  "ident_percent_std": sample_std_ident_percent}
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


def process_info_sample(ident_details, experiment_details, experiment_type_indices,
                        ident_fun_choice=(), combine_fluxes=0):
    if ident_fun_choice:
        if len(ident_fun_choice) <= 1:
            combine_fluxes = 0

    print("Process information From Identifiability Analysis.....\n")
    number_of_samples = len(ident_details)
    all_sample_data_list = []
    all_sample_max_parameter = []
    all_sample_true_parameter = []
    all_sample_experiment_list = []
    all_sample_combined_flux_data_list = []
    all_sample_combined_flux_max_parameter = []
    all_sample_combined_flux_true_parameter = []
    all_sample_combined_flux_experiment_info = []
    number_of_fluxes_per_sample = []
    for j_sample, j_sample_ident_detail in enumerate(ident_details):
        print("Processing identifiability data for sample {} of {}".format(j_sample+1, number_of_samples))
        # collect flux based identifiability information/data
        all_flux_data_list, all_flux_max_parameter, \
            all_flux_true_parameter = flux_based_ident_info(j_sample_ident_detail)
        all_sample_data_list.append(all_flux_data_list)
        all_sample_max_parameter.append(all_flux_max_parameter)
        all_sample_true_parameter.append(all_flux_true_parameter)
        number_of_fluxes_per_sample.append(len(all_flux_data_list))

        # collect flux based information on experiments contributing to identifiability
        all_flux_experiment_list = flux_based_experiment_info(j_sample_ident_detail,
                                                              experiment_details[j_sample],
                                                              experiment_type_indices)
        all_sample_experiment_list.append(all_flux_experiment_list)

        # collate data for all fluxes using the same combination of experimental data
        if combine_fluxes:
            print("Processing identifiability data combining all fluxes in sample")
            # collate all identifiability information
            combined_flux_ident_data = collate_flux_based_data(j_sample_ident_detail)
            # process/collect corresponding data for all fluxes using the same combination of experimental data
            combined_flux_data_list, combined_flux_max_parameter, \
                combined_flux_true_parameter = process_info(combined_flux_ident_data)
            all_sample_combined_flux_data_list.append(combined_flux_data_list)
            all_sample_combined_flux_max_parameter.append(combined_flux_max_parameter)
            all_sample_combined_flux_true_parameter.append(combined_flux_true_parameter)
            # collect experiment information for all fluxes using the same combination of experimental data
            combined_flux_experiment_info = experiments_in_ident_data(combined_flux_ident_data["boolean"],
                                                                      experiment_details[j_sample],
                                                                      experiment_type_indices,
                                                                      combined_flux_ident_data["flux id"],
                                                                      combined_flux_ident_data["flux choice"])
            all_sample_combined_flux_experiment_info.append(combined_flux_experiment_info)
            print("Combined flux identifiability analysis complete \n")

    # generate averages and standard deviations for multi sample data combinations for individual fluxes
    processed_all_flux_data_list, \
        processed_all_flux_max_parameter, \
        processed_all_flux_experiment_info, \
        processed_all_flux_parameter_values = sample_based_averages(number_of_fluxes_per_sample, all_sample_data_list,
                                                                    all_sample_max_parameter,
                                                                    all_sample_experiment_list,
                                                                    all_sample_true_parameter)
    all_sample_data_utility = {"raw": all_sample_data_list,
                               "processed": processed_all_flux_data_list}
    all_sample_parameter_identifiability = {"raw": all_sample_max_parameter,
                                            "processed": processed_all_flux_max_parameter}
    all_sample_parameter_value = {"raw": all_sample_true_parameter,
                                  "processed": processed_all_flux_parameter_values}
    all_sample_experiment_info = {"raw": all_sample_experiment_list,
                                  "processed": processed_all_flux_experiment_info}
    if combine_fluxes:
        processed_all_sample_combined_flux_data_list = \
            combined_sample_based_averages_data_utility(all_sample_combined_flux_data_list)
        all_sample_combined_data_utility = {"raw": all_sample_combined_flux_data_list,
                                            "processed": processed_all_sample_combined_flux_data_list}
        all_sample_combined_parameter_identifibaility = {"raw": all_sample_combined_flux_max_parameter,
                                                         "processed": []}
        processed_all_sample_combined_parameter_value = \
            combined_sampled_based_averages_parameter_value(all_sample_combined_flux_true_parameter)
        all_sample_combined_parameter_value = {"raw": all_sample_combined_flux_true_parameter,
                                               "processed": processed_all_sample_combined_parameter_value}
        processed_all_sample_combined_flux_experiment_info = \
            combined_sample_based_averages_experiment_info(all_sample_combined_flux_experiment_info)
        all_sample_combined_experiment_info = {"raw": all_sample_combined_flux_experiment_info,
                                               "processed": processed_all_sample_combined_flux_experiment_info}

        return all_sample_data_utility, all_sample_parameter_identifiability, all_sample_parameter_value, \
               all_sample_experiment_info, all_sample_combined_data_utility, \
               all_sample_combined_parameter_identifibaility, all_sample_combined_parameter_value, \
               all_sample_combined_experiment_info
    else:
        return all_sample_data_utility, all_sample_parameter_identifiability, \
               all_sample_parameter_value, all_sample_experiment_info, [], [], [], []


def extract_parameter_values(parameter_value):
    """extract all parameter values in a given flux to re-simulate model with newly determined parameters"""
    processed_info = parameter_value["processed"]
    number_flux = len(processed_info)
    try:
        for i_flux, i_flux_info in enumerate(processed_info):
            print('\nExtracting parameters for flux {} of {}'.format(i_flux + 1, number_flux))
            number_samples, number_data = i_flux_info[0]["raw sample ident data"].shape
            # collect parameter names
            all_parameter_name = [i_parameter_info["parameter name"] for i_parameter_info in i_flux_info]
            # collect boolean info for all parameters
            all_parameter_boolean_info = [i_parameter_info["raw sample boolean data"]
                                          for i_parameter_info in i_flux_info]
            # collect parameter value info for all parameters
            all_parameter_ident_info = [i_parameter_info["raw sample ident data"] for i_parameter_info in i_flux_info]
            # collect parameter values from each sample
            all_sample_info = []
            all_sample_boolean = []
            for j_sample in range(0, number_samples):
                j_sample_info = [tuple([np.array(i_parameter_ident[j_sample, i_data])
                                        for i_parameter_ident in all_parameter_ident_info])
                                 for i_data in range(0, number_data)]
                j_sample_boolean = [tuple([i_parameter_ident[j_sample][i_data]
                                           for i_parameter_ident in all_parameter_boolean_info])
                                    for i_data in range(0, number_data)]
                all_sample_info.append(j_sample_info)
                all_sample_boolean.append(j_sample_boolean)
            # get parameters values only when all parameters are identified
            all_sample_ident_info = []
            for k_sample, (k_sample_boolean, k_sample_info) in enumerate(zip(all_sample_boolean, all_sample_info)):
                boolean_sum = [all(j_data_boolean) for j_data_boolean in k_sample_boolean]
                all_true_ident = [(j_data_id, j_data_info)
                                  for j_data_id, j_data_info in enumerate(k_sample_info) if boolean_sum[j_data_id]]
                k_sample_parameter_value_info = {"data_id": [data_id_value[0] for data_id_value in all_true_ident],
                                                 "parameter_value": [dict(zip(all_parameter_name, data_id_value[1]))
                                                                     for data_id_value in all_true_ident],
                                                 "parameter name": all_parameter_name}
                all_sample_ident_info.append(k_sample_parameter_value_info)
    except TypeError:
        all_sample_ident_info = []
    return all_sample_ident_info


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
