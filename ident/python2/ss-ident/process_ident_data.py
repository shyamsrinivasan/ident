import numpy as np
import itertools as it


def parameter_based_processing(ident_details):
    """parameter-based classification of experimental datasets"""
    number_data, p = ident_details["boolean"].shape
    # p = np.cumsum([0] + parameters_per_flux).tolist()
    fp_list_keys = ['flux{}p{}'.format(f_index + 1, p_index + 1)
                    for f_index, p_limit in enumerate(ident_details["parameters"])
                    for p_index in range(0, p_limit)]
    fp_list_values = [[data_id for data_id in it.compress(range(0, number_data), list(ident_details["boolean"][:, parameter_id]))]
                      for parameter_id in range(0, p)]
    fp_list = dict(zip(fp_list_keys, fp_list_values))
    return fp_list


def data_based_processing(ident_details):
    """dataset dependent classification of parameters"""
    number_data, p = ident_details["boolean"].shape
    data_list_keys = ['dataset{}'.format(perturbation_id + 1) for perturbation_id in range(0, number_data)]
    data_list_values = [[parameter_id for parameter_id in it.compress(range(0, p), list(ident_details["boolean"][data_id, :]))]
                        for data_id in range(0, number_data)]
    data_list = dict(zip(data_list_keys, data_list_values))
    return data_list


def process_ident_data(ident_values, number_data):
    # create signed boolean array for identifiability
    signed_ident_values = np.sign(ident_values)
    ident_fun_val = []
    for id in range(0, number_data):
        ident_fun_val.append(signed_ident_values[id * 12:(id + 1) * 12, -1])
    p_list = [[p_id for p_id, val in enumerate(data_set) if val > 0] for data_set in ident_fun_val]
    p_list_boolean = [[True if parameter_id in list_1 else False for parameter_id in range(0, 12)] for list_1 in p_list]
    return p_list, np.array(p_list_boolean)


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


def get_useful_data_info(ident_boolean_array, experiment_details, perturbation_details, data_needed):
    """get info all data sets that identify maximum number of parameters"""
    _, p = ident_boolean_array.shape
    # data_needed = max_data["id"]
    parameters_identified_by_max_data_id = []
    # get parameters identified by data_needed
    data_needed_id = 0
    for data_id in data_needed:
        # data_set_name = 'option {}'.format(data_needed_id + 1)
        # data_set_ident = [ident_parameter_name(i)
        #                   for i in it.compress(range(0, p), list(ident_boolean_array[data_id, :]))]
        data_set_ident_id = [i for i in it.compress(range(0, p), list(ident_boolean_array[data_id, :]))]
        # identify all experiments for given data set
        #perturbation_id = [int(i) for i in experiment_details["details"][data_id, range(0, 12, 4)].tolist()]
        perturbation_id = experiment_details["experiment_id"][data_id]
        # identify parameters perturbed by said experiment
        # perturbed_parameters = perturbation_details["indices"][perturbation_id, :]
        # perturbed_parameters = experiment_details["parameter_ids"][data_id]
        # get parameter names from parameter indices
        # perturbed_parameter_name = [(model_parameter_name(int(j[0])), j[1]) for j in perturbed_parameters.tolist()]
        info = {'id': data_id, 'parameter_ids':data_set_ident_id,
                'experiment_id':perturbation_id, 'ssid':experiment_details["ssid"][data_id]}
        parameters_identified_by_max_data_id.append(info)
        data_needed_id += 1
    #data_list = dict(parameters_identified_by_max_data_id)
    return parameters_identified_by_max_data_id


def get_additional_data(ident_details, experiment_details, perturbation_details, unidentified_parameter_ids):
    """find data sets that identify parameters in unidentified set"""
    number_data, p = ident_details["boolean"].shape
    # get boolean for unidentified parameters
    unident_boolean = ident_details["boolean"][:, unidentified_parameter_ids]
    # work on this later
    # max_data = get_most_useful_dataset(unident_boolean)
    # if max_data["maximum"] > 1:
    #    selected_data_id = max_data["id"]
    #else:

    # get every single data set that identifies unidentified parameters
    data_ident = [[k for k in it.compress(range(0, number_data), list(unident_boolean[:, unid_p]))]
                  for unid_p in range(0, len(unidentified_parameter_ids))]
    # get other info for data sets in this list
    extra_data_details = []
    for list_elements in data_ident:
        if list_elements:
            data_list = get_useful_data_info(ident_details["boolean"],
                                             experiment_details,
                                             perturbation_details,
                                             list_elements)
            for dictionary in data_list:
                extra_data_details.append(dictionary)
        else:
            extra_data_details.append({})



    # collect common data ids in data_ident
    #for j_data_ident in data_ident:
    #for i_data in range(0, number_data):
    #   unidentifiable_parameters = [k for k in it.compress(range(0, p), list(np.logical_not(ident_details["boolean"][i_data, :])))]

    #    identified_parameters = [i for i in it.compress(range(0, p), list(ident_details["boolean"][i_data, :]))]
    #    identifiable_parameter_boolean = [True if un_id in identified_parameters
    #                                      else False
    #                                      for un_id in unidentified_parameter_ids]
    #    if any(identifiable_parameter_boolean):
    #        identifiable_parameter = [j
    #                                  for j in it.compress(unidentified_parameter_ids, identifiable_parameter_boolean)]
    #        extra_data_details.append({'data':i_data, 'parameter':identifiable_parameter})
    return extra_data_details


def calculate_experiment_combos(ident_details, experiment_details, perturbation_details, data_list):
    """choose additional data sets and consequently, experiments (if possible) to identify other parameters not
    identified by chosen data set(s)"""
    # get parameters identified by chosen data set
    number_data, p = ident_details["boolean"].shape
    new_combos = []
    new_combo_number = 0
    for chosen_data in data_list:
        # chosen_data = data_list[option_id]
        data_set_id = chosen_data["id"]
        # get parameters identified
        identified_parameters = [i for i in it.compress(range(0, p), list(ident_details["boolean"][data_set_id, :]))]
        # get unidentified parameters
        nonidentified_parameters = list(set(range(0, p)) - set(identified_parameters))

        # go through all data sets to identify indexes that identify unidentified parameters
        extra_data = get_additional_data(ident_details,
                                         experiment_details,
                                         perturbation_details,
                                         nonidentified_parameters)
        for extra_id, dict_id in enumerate(extra_data):
            if dict_id:
                new_option_name = 'combination {}'.format(new_combo_number + 1)
                new_option_composition = [chosen_data["id"], dict_id["id"]]
                new_option_experiment_id = [chosen_data["experiment_id"], dict_id["experiment_id"]]
                new_option_ssid = [chosen_data["ssid"], dict_id["ssid"]]
                new_options_parameter_ids = \
                    list(set(chosen_data["parameter_ids"]) | set(dict_id["parameter_ids"]))
                combo_info = {'id': new_option_composition, 'parameter_ids': new_options_parameter_ids,
                              'experiment_ids': new_option_experiment_id,
                              'ssid': new_option_ssid}
                new_combos.append(combo_info)
                new_combo_number += 1
    return new_combos


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
    percentage_of_total = []
    for i, j_parameter_identified in enumerate(number_parameters_identified):
        number_data_identifying_j = len(data_id[i])
        percentage_of_total.append(float(number_data_identifying_j) / float(number_data) * 100)
    data_usefulness = {'number_parameters_ided': number_parameters_identified,
                       'index': data_id,
                       'percentage': percentage_of_total}
    return data_usefulness


def process_info(ident_details, experiment_details, perturbation_details):
    number_data, p = ident_details["boolean"].shape

    # get data identification percentages to classify utility of data sets
    data_usefulness = data_usefulness_percentage(ident_details)

    # get info all the aforementioned data sets
    original_data = []
    combination_data = []
    for i_data in range(0, len(data_usefulness["number_parameters_ided"])):
        # get all info on all data sets present in useful data from above
        temp_list = get_useful_data_info(ident_details["boolean"],
                                         experiment_details,
                                         perturbation_details,
                                         data_usefulness["index"][i_data])
        original_data.append(temp_list)

        # choose additional data sets and consequently, experiments (if possible) to identify other parameters not
        # identified by chosen data set(s)
        new_combos = calculate_experiment_combos(ident_details,
                                                 experiment_details,
                                                 perturbation_details,
                                                 temp_list)
        combination_data.append(new_combos)

    # print total individual data sets required and combinations found
    number_original_data = 0
    number_combination_data = 0
    for list_pos, i_data in enumerate(original_data):
        # original data sets
        print('Original Data sets that can detect {} parameters: {}'.
              format(data_usefulness["number_parameters_ided"][list_pos], len(i_data)))
        number_original_data += len(i_data)
        # combination data sets
        print('Combination Data sets for more parameters: {}'.format(len(combination_data[list_pos])))
        number_combination_data += len(combination_data[list_pos])

    # decide which experiments to perform for each parameter based on above calculations

    # most easily identifieable parameter - based on frequency of identification
    identifying_data = [] # np.zeros((number_parameter, 1))
    for iparameter in range(0, p):
        idata = [data_id for data_id in it.compress(range(0, number_data), ident_details["boolean"][:, iparameter])]
        identifying_data.append(len(idata))

    # get parameters identified by most data sets (most identifiable parameter)
    max_identified = max(identifying_data)
    max_parameter_id = [identifying_data.index(max(identifying_data))]
    for iparameter, number_identifying_data in enumerate(identifying_data):
        if number_identifying_data == max_identified and iparameter not in max_parameter_id:
            max_parameter_id.append(iparameter)
        elif number_identifying_data > max_identified:
            max_identified = number_identifying_data
            max_parameter_id = [iparameter]
    max_parameter = {"maximum": max_identified,
                     "id": max_parameter_id,
                     "data": identifying_data}

    # ident_parameter_names = ident_parameter_name(range(0, 12))

    return data_usefulness, original_data, combination_data, max_parameter
