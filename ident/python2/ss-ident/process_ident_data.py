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
    max_data_id = [identified_parameters.index(max(identified_parameters))]
    for idata, number_identified_parameter in enumerate(identified_parameters):
        if number_identified_parameter == max_useful_data and idata not in max_data_id:
            max_data_id.append(idata)
        elif number_identified_parameter > max_useful_data:
            max_useful_data = number_identified_parameter
            max_data_id = [idata]
    max_data = {"maximum": max_useful_data,
                "id": max_data_id,
                "parameters": identified_parameters}
    return max_data


def get_useful_data_info(ident_boolean_array, experiment_details, perturbation_details, data_needed,
                         ident_parameter_name, model_parameter_name):
    """get info all data sets that identify maximum number of parameters"""
    _, p = ident_boolean_array.shape
    # data_needed = max_data["id"]
    parameters_identified_by_max_data_id = []
    # get parameters identified by data_needed
    data_needed_id = 0
    for data_id in data_needed:
        data_set_name = 'option {}'.format(data_needed_id + 1)
        data_set_ident = [ident_parameter_name(i)
                          for i in it.compress(range(0, p), list(ident_boolean_array[data_id, :]))]
        data_set_ident_id = [i for i in it.compress(range(0, p), list(ident_boolean_array[data_id, :]))]
        # identify all experiments for given data set
        #perturbation_id = [int(i) for i in experiment_details["details"][data_id, range(0, 12, 4)].tolist()]
        perturbation_id = experiment_details["experiment_id"][data_id]
        # identify parameters perturbed by said experiment
        perturbed_parameters = perturbation_details["indices"][perturbation_id, :]
        #perturbed_parameters = experiment_details["parameter_ids"][data_id]
        # get parameter names from parameter indices
        perturbed_parameter_name = [(model_parameter_name(int(j[0])), j[1]) for j in perturbed_parameters.tolist()]
        info = {'id': data_id, 'parameters': data_set_ident,
                'experiments': perturbed_parameter_name, 'parameter_ids':data_set_ident_id,
                'experiment_id':perturbation_id, 'ssid':experiment_details["ssid"][data_id]}
        parameters_identified_by_max_data_id.append((data_set_name, info))
        data_needed_id += 1
    data_list = dict(parameters_identified_by_max_data_id)
    return data_list


def get_additional_data(ident_details, experiment_details, perturbation_details,
                        unidentified_parameter_ids, ident_parameter_name, model_parameter_name):
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
                                             list_elements,
                                             ident_parameter_name, model_parameter_name)
            extra_data_details.append(data_list)
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


def calculate_experiment_combos(ident_details, experiment_details, perturbation_details, data_list,
                                ident_parameter_name, model_parameter_name):
    """choose additional data sets and consequently, experiments (if possible) to identify other parameters not
    identified by chosen data set(s)"""
    # get parameters identified by chosen data set
    number_data, p = ident_details["boolean"].shape
    new_combos = []
    new_combo_number = 0
    for option_id in data_list:
        chosen_data = data_list[option_id]
        data_set_id = chosen_data["id"]
        # get parameters identified
        identified_parameters = [i for i in it.compress(range(0, p), list(ident_details["boolean"][data_set_id, :]))]
        # get unidentified parameters
        nonidentified_parameters = list(set(range(0, p)) - set(identified_parameters))

        # go through all data sets to identify indexes that identify unidentified parameters
        extra_data = get_additional_data(ident_details,
                                         experiment_details,
                                         perturbation_details,
                                         nonidentified_parameters,
                                         ident_parameter_name,
                                         model_parameter_name)
        for extra_id, dict_id in enumerate(extra_data):
            for id, options in enumerate(dict_id):
                new_option_name = 'combination {}'.format(new_combo_number + 1)
                new_option_composition = [chosen_data["id"], dict_id[options]["id"]]
                new_option_experiments = [chosen_data["experiments"], dict_id[options]["experiments"]]
                new_options_parameter_ids = list(
                    set(chosen_data["parameter_ids"]) | set(dict_id[options]["parameter_ids"]))
                new_option_parameters = ident_parameter_name(new_options_parameter_ids)
                combo_info = {'id': new_option_composition,
                              'experiments': new_option_experiments, 'parameters': new_option_parameters,
                              'parameter_ids': new_options_parameter_ids}
                new_combos.append((new_option_name, combo_info))
                new_combo_number += 1
    return new_combos


def process_info(ident_details, experiment_details, perturbation_details, number_fluxes,
                 ident_parameter_name, model_parameter_name):
    number_data, p = ident_details["boolean"].shape
    # parameter-based classification of experimental datasets
    fp_list = parameter_based_processing(ident_details)

    # dataset dependent classification of parameters
    # data_list = data_based_processing(ident_details)

    # most useful dataset - based on number of parameter identified
    max_data = get_most_useful_dataset(ident_details["boolean"])

    # decide which experiments to perform for each parameter based on above calculations
    # get all info on most identifiable data sets (most useful data set)
    data_list = get_useful_data_info(ident_details["boolean"],
                                     experiment_details,
                                     perturbation_details,
                                     max_data["id"],
                                     ident_parameter_name, model_parameter_name)

    # choose additional data sets and consequently, experiments (if possible) to identify other parameters not
    # identified by chosen data set(s)
    new_combos = calculate_experiment_combos(ident_details, experiment_details, perturbation_details, data_list,
                                             ident_parameter_name, model_parameter_name)

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

    return data_list, new_combos, max_parameter
