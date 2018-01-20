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


def get_useful_data_info(ident_boolean_array, experiment_details, perturbation_details, data_needed):
    """get info all data sets that identify maximum number of parameters"""
    _, p = ident_boolean_array.shape
    parameters_identified_by_max_data_id = []
    # get parameters identified by data_needed
    all_ided_parameters = parameters_ided_by_combination(ident_boolean_array, data_needed)
    # get experiment ids used in/that form data combination
    all_experiment_ids = experiments_in_combination(data_needed, experiment_details)

    data_needed_id = 0
    for j_id, data_id in enumerate(data_needed):
        data_set_ident_id = all_ided_parameters[j_id]
        # identify all experiments for given data set
        #perturbation_id = [int(i) for i in experiment_details["details"][data_id, range(0, 12, 4)].tolist()]
        perturbation_id = all_experiment_ids[j_id]
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
    total_number_of_data = []
    percentage_of_total = []
    for i, j_parameter_identified in enumerate(number_parameters_identified):
        number_data_identifying_j = len(data_id[i])
        total_number_of_data.append(number_data_identifying_j)
        percentage_of_total.append(float(number_data_identifying_j) / float(number_data) * 100)
    data_usefulness = {'number': number_parameters_identified,  # number of parameters ided
                       'index': data_id,                        # index of data combinations that id x parameters
                       'total': total_number_of_data,           # number of data combinations that id x parameters
                       'percentage': percentage_of_total,       # percentage of data combinations that id x parameters
                       'number of data combinations': number_data,  # total number of data combinations used
                       'flux id': ident_details["flux id"]}     # flux id of data being processed
    return data_usefulness


def flux_parameter_plot_data(original_data, number_of_parameters, case=1):
    all_boolean_p_id = []
    for len_pos, i_list in enumerate(original_data):
        for i_data in i_list:
            boolean_p_id = [True if j_p in i_data["parameter_ids"] else False for j_p in range(0, number_of_parameters)]
            all_boolean_p_id.append(boolean_p_id)
    number_data, number_p = np.array(all_boolean_p_id).shape
    all_boolean_p_id = [list(j_p) for j_p in np.transpose(np.array(all_boolean_p_id))]
    if case == 1:
        # get total data identifying each parameter
        all_boolean_p_id_sum = [sum(j_list) for j_list in all_boolean_p_id]
        all_boolean_p_id_fraction = [float(sum(j_list))*100/number_data for j_list in all_boolean_p_id]
        return all_boolean_p_id_sum, all_boolean_p_id_fraction, all_boolean_p_id
    elif case == 2:
        all_boolean_e_id = []
        for len_pos, i_list in enumerate(original_data):
            for i_data in i_list:
                # experiments done
                boolean_e_id = [True if j_p in i_data["experiment_id"] else False for j_p in range(0, 18)]
                all_boolean_e_id.append(boolean_e_id)
        # all_boolean_e_id = [list(j_p) for j_p in np.transpose(np.array(all_boolean_e_id))]
        # get total for each experiment in each identifying data set
        # all_boolean_e_id = [sum(k_list) for k_list in all_boolean_e_id]
        return all_boolean_p_id, all_boolean_e_id
    else:
        return []


def single_flux_parameter_data(ident_data, number_of_parameters_per_flux, case):
    number_of_fluxes = len(ident_data)
    # all_flux_total_p = []
    all_flux_info = []
    for iflux, iflux_data in enumerate(ident_data):
        total_p, fraction_p, all_p_boolean = flux_parameter_plot_data(iflux_data,
                                                                      number_of_parameters_per_flux[iflux],
                                                                      case)
        all_flux_info.append({"total": total_p, "fraction": fraction_p, "boolean": all_p_boolean})
    return all_flux_info


def parameter_plot_data_per_sample(original_data, number_of_parameters_per_flux, case=1):
    number_of_samples = len(original_data)
    all_sample_total_p = []
    all_sample_fraction_p = []
    all_sample_all_p_boolean = []
    for i_sample_ident_data in original_data:
        all_flux_info = single_flux_parameter_data(i_sample_ident_data, number_of_parameters_per_flux, case)
        # total_p, fraction_p, all_p_boolean = flux_parameter_plot_data(i_sample_ident_data, case)
        all_sample_total_p.append(total_p)
        all_sample_fraction_p.append(fraction_p)
        all_sample_all_p_boolean.append(all_p_boolean)

    # collect data from all samples and calculate means and standard deviations for each parameter
    all_sample_total_means = list(np.mean(np.array(all_sample_total_p), axis=0))
    all_sample_total_std = list(np.std(np.array(all_sample_total_p), axis=0))
    all_sample_fraction_means = list(np.mean(np.array(all_sample_fraction_p), axis=0))
    all_sample_fraction_std = list(np.std(np.array(all_sample_fraction_p), axis=0))
    all_sample_totals = {'means': all_sample_total_means, 'std': all_sample_total_std}
    all_sample_fractions = {'means': all_sample_fraction_means, 'std': all_sample_fraction_std}
    for j_sample in range(0, number_of_samples):
        pass
    return all_sample_totals, all_sample_fractions, all_sample_all_p_boolean


def experiments_in_ident_data(boolean_ident_data, experiment_data, experiment_type_index, flux_id):
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
            experiment_type_boolean = [[True if exp_id in i_experiment_type else False
                                        for exp_id in experiments_in_identifying_data[:, j_position_in_combination]]
                                       for i_experiment_type in experiment_type_index]

            # get experiment type frequency
            experiment_type_frequency = [sum(i_experiment_type_boolean)
                                         for i_experiment_type_boolean in experiment_type_boolean]
            total_identifying_data = sum(experiment_type_frequency)
            experiment_type_percentage = [float(i_frequency)*100/total_identifying_data
                                          for i_frequency in experiment_type_frequency]

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
                                                 "flux id": flux_id})
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
                                                                   flux_id=j_flux_detail["flux id"])
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
                     "flux id": ident_details["flux id"]}
    return max_parameter


def process_info(ident_details):
    """get data utility and parameter identifiability and other data for
    parameters of each flux passed as input"""
    # get data identification percentages to classify utility of data sets
    data_list = data_utility(ident_details)
    # most easily identifiable parameter - based on frequency of identification
    max_parameter = parameter_identifiability(ident_details)

    return data_list, max_parameter


def flux_based_ident_info(sample_ident_detail):
    """parse information by looping though each flux present within each sample that is passed an input argument"""
    number_of_fluxes = len(sample_ident_detail)
    all_flux_data_list = []
    all_flux_max_parameter = []
    for j_flux, j_flux_info in enumerate(sample_ident_detail):
        print("Processing identifiability for flux {} of {}".format(j_flux + 1, number_of_fluxes))
        data_list, max_parameter = process_info(j_flux_info)
        all_flux_data_list.append(data_list)
        all_flux_max_parameter.append(max_parameter)
        print("Information Processing Complete for flux {} \n".format(j_flux + 1))
    return all_flux_data_list, all_flux_max_parameter


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
    combined_flux_ident_data = {"boolean": all_flux_boolean,
                                "values": all_flux_values,
                                "parameters": all_flux_parameters,
                                "flux id": flux_id_list}
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
                                  "flux id": all_sample_data_list[0][j_flux]["flux id"]}
        j_flux_processed_percent = {"mean": j_flux_percent_mean,
                                    "std": j_flux_percent_std,
                                    "number": j_flux_number_mean,
                                    "flux id": all_sample_data_list[0][j_flux]["flux id"]}
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
                                  "flux id": all_sample_max_parameter[0][j_flux]["flux id"]}
        j_flux_processed_percent = {"mean": j_flux_percent_mean,
                                    "std": j_flux_percent_std,
                                    "flux id": all_sample_max_parameter[0][j_flux]["flux id"]}
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
                                             [k_parameter][i_position]["flux id"]}
                i_position_percent_info = {"mean": i_position_percent_mean,
                                           "std": i_position_percent_std,
                                           "flux id": all_sample_experiment_list[0][j_flux]
                                           [k_parameter][i_position]["flux id"]}
                # collect data on each position
                all_position_info.append({"total": i_position_frequency_info,
                                          "percentage": i_position_percent_info,
                                          "position": i_position})
            # collect info on each parameter in all positions
            all_parameter_all_position_info.append(all_position_info)
        # collect data on all parameters for each flux
        all_flux_experiment_info.append(all_parameter_all_position_info)

    return all_flux_experiment_info


def sample_based_averages(number_of_fluxes_per_sample,
                          all_sample_data_list, all_sample_max_parameter, all_sample_experiment_list):
    """call sample based collate functions and generate averages and standard deviations for
    both data utility and parameter identifiability for each flux"""
    all_flux_data_utility = collate_sample_based_data_utility(number_of_fluxes_per_sample,
                                                              all_sample_data_list)
    all_flux_max_parameter = collate_sample_based_identifibaility(number_of_fluxes_per_sample,
                                                                  all_sample_max_parameter)
    all_flux_experiment_info = collate_sample_based_experiment_info(number_of_fluxes_per_sample,
                                                                    all_sample_experiment_list)
    return all_flux_data_utility, all_flux_max_parameter, all_flux_experiment_info


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
                           "flux id": all_sample_combined_flux_data_list[0]["flux id"]}
    all_processed_percent = {"mean": sample_percent_mean,
                             "std": sample_percent_std,
                             "number": sample_number_mean,
                             "total data": all_sample_combined_flux_data_list[0]["number of data combinations"],
                             "flux id": all_sample_combined_flux_data_list[0]["flux id"]}
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
                                         [j_parameter][i_position]["flux id"]}
            i_position_percent_info = {"mean": i_position_percent_mean,
                                       "std": i_position_percent_std,
                                       "flux id": all_sample_combined_flux_experiment_info[0]
                                       [j_parameter][i_position]["flux id"]}
            # collect data on each position
            all_position_info.append({"frequency": i_position_frequency_info,
                                      "percentage": i_position_percent_info,
                                      "position": i_position})
        # collect data on each parameter
        all_parameter_all_position_info.append(all_position_info)

    return all_parameter_all_position_info


def process_info_sample(ident_details, experiment_details, experiment_type_indices, perturbation_details, combine_fluxes=0):
    print("Process information From Identifiability Analysis.....\n")
    number_of_samples = len(ident_details)
    all_sample_data_list = []
    all_sample_max_parameter = []
    all_sample_experiment_list = []
    all_sample_combined_flux_data_list = []
    all_sample_combined_flux_max_parameter = []
    all_sample_combined_flux_experiment_info = []
    number_of_fluxes_per_sample = []
    for j_sample, j_sample_ident_detail in enumerate(ident_details):
        print("Processing identifiability data for sample {} of {}".format(j_sample+1, number_of_samples))
        # collect flux based identifiability information/data
        all_flux_data_list, all_flux_max_parameter = flux_based_ident_info(j_sample_ident_detail)
        all_sample_data_list.append(all_flux_data_list)
        all_sample_max_parameter.append(all_flux_max_parameter)
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
            combined_flux_data_list, combined_flux_max_parameter = process_info(combined_flux_ident_data)
            all_sample_combined_flux_data_list.append(combined_flux_data_list)
            all_sample_combined_flux_max_parameter.append(combined_flux_max_parameter)
            # collect experiment information for all fluxes using the same combination of experimental data
            combined_flux_experiment_info = experiments_in_ident_data(combined_flux_ident_data["boolean"],
                                                                      experiment_details[j_sample],
                                                                      experiment_type_indices,
                                                                      combined_flux_ident_data["flux id"])
            all_sample_combined_flux_experiment_info.append(combined_flux_experiment_info)
            print("Combined flux identifiability analysis complete \n")

    # generate averages and standard deviations for multi sample data combinations for individual fluxes
    processed_all_flux_data_list, \
    processed_all_flux_max_parameter,\
    processed_all_flux_experiment_info = sample_based_averages(number_of_fluxes_per_sample,
                                                               all_sample_data_list,
                                                               all_sample_max_parameter,
                                                               all_sample_experiment_list)
    all_sample_data_utility = {"raw": all_sample_data_list,
                               "processed": processed_all_flux_data_list}
    all_sample_parameter_identifiability = {"raw": all_sample_max_parameter,
                                            "processed": processed_all_flux_max_parameter}
    all_sample_experiment_info = {"raw": all_sample_experiment_list,
                                  "processed": processed_all_flux_experiment_info}
    if combine_fluxes:
        processed_all_sample_combined_flux_data_list = \
            combined_sample_based_averages_data_utility(all_sample_combined_flux_data_list)
        all_sample_combined_data_utility = {"raw": all_sample_combined_flux_data_list,
                                            "processed": processed_all_sample_combined_flux_data_list}
        all_sample_combined_parameter_identifibaility = {"raw": all_sample_combined_flux_max_parameter,
                                                         "processed": []}
        processed_all_sample_combined_flux_experiment_info = \
            combined_sample_based_averages_experiment_info(all_sample_combined_flux_experiment_info)
        all_sample_combined_experiment_info = {"raw": all_sample_combined_flux_experiment_info,
                                               "processed": processed_all_sample_combined_flux_experiment_info}
        return all_sample_data_utility, all_sample_parameter_identifiability, all_sample_experiment_info, \
               all_sample_combined_data_utility, all_sample_combined_parameter_identifibaility, \
               all_sample_combined_experiment_info
    else:
        return all_sample_data_utility, all_sample_parameter_identifiability, all_sample_experiment_info


def get_data_combinations(original_data, chosen_data_id):
    """get new data combinations that result in additional parameters being identified
    for a given data set until all parameters can be identified"""
    # given - chosen_data_id - id of data sets already chosen
    # return - additional data sets/experiments to enable identification of parameters not identified by chosen_data_id
    temp_list = get_useful_data_info(ident_details["boolean"],
                                     experiment_details,
                                     perturbation_details,
                                     data_usefulness["index"][i_data])
    new_combos = calculate_experiment_combos(ident_details,
                                             experiment_details,
                                             perturbation_details,
                                             temp_list)
    return None
