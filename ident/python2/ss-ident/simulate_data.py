import numpy as np
import itertools as it


def parameter_change(new_value, old_value):
    return (new_value - old_value) *100 / old_value


def get_changed_parameters(original_parameters, changed_parameters, experiment_index, parameter_index):
    """create dictionary of parameter values after perturbations -
    determine parameter changes and percentage changes"""
    parameter_list = {"experiment index": experiment_index,
                      "parameter name": parameter_index,
                      "parameter value": changed_parameters[experiment_index][parameter_index],
                      "parameter change": parameter_change(changed_parameters[experiment_index][parameter_index],
                                                           original_parameters[parameter_index])}
    return parameter_list


def get_data_combinations(xss, experiments_per_set, experiment_choice, combination_choice):
    """get combinations of data sets based on all perturbations performed on the model"""

    # select only experiment indices provided in experiment_choice
    total_number_of_experiments = len(xss)
    if experiment_choice:
        chosen_xss = [xss[i_experiment] for i_experiment in experiment_choice]
        number_of_experiments = len(chosen_xss)
    else:
        chosen_xss = xss
        number_of_experiments = len(chosen_xss)
        experiment_choice = range(0, number_of_experiments)
    data_combinations = []
    # get all experiment combinations based on experiments_per_set
    for item in it.permutations(experiment_choice, experiments_per_set):
        data_combinations.append(item)

    # choose only combination_choice values of experimental data set combinations
    if combination_choice:
        try:
            data_combinations = data_combinations[:combination_choice]
            combination_choice = range(0, combination_choice)
        except TypeError:
            new_data_combinations = []
            for indexes in combination_choice:
                new_data_combinations.append(data_combinations[indexes])
            data_combinations = new_data_combinations[:]
    else:
        number_of_combinations = len(data_combinations)
        combination_choice = range(0, number_of_combinations)

    # create boolean numpy array of perturbations in each data set
    data_combination_boolean = [[True if experiment_id in list_1 else False
                                 for experiment_id in range(0, total_number_of_experiments)]
                                for list_1 in data_combinations]
    data_combination_boolean = np.array(data_combination_boolean)
    return data_combinations, data_combination_boolean, experiment_choice, combination_choice


def data_for_each_sample(perturbation_details, experiments_per_set,
                         data_combinations, xss, fss, flux_id, choose):
    """get simulated experimental data for each noisy sample supplied as input argument"""
    # collect data for combinations in data_combinations for each sample in xss
    parameters = perturbation_details["values"]
    experiment_indices = perturbation_details["indices"]
    original = perturbation_details["original"]
    ssid = perturbation_details["ssid"]

    # initialize vectors/arrays and other lists for arrangement
    number_data_sets = len(data_combinations)
    size_of_data_set = 8
    dataset_value_array = np.zeros((number_data_sets, size_of_data_set * experiments_per_set))
    # all_parameter_change = np.zeros((number_data_sets, 4 * experiments_per_set))
    all_parameter_change = []

    # arrange and collect parameter/ss values
    experiment_index_iter = []
    all_parameter_ids = []
    all_parameter_values = []
    all_parameter_changes = []
    all_ssid = []
    for choice_index, data_index in zip(choose, enumerate(data_combinations)):
        data = []
        parameter_id = []
        changes = {}
        experiment_index_iter.append(data_index[1])
        parameter_value = []
        parameter_change = []
        ssid_changes = []
        for iter, index in enumerate(data_index[1]):
            data.append(np.hstack((parameters[index]["ac"],
                                   xss[index],
                                   fss[index][flux_id])))
            parameter_list = get_changed_parameters(original_parameters=original,
                                                    changed_parameters=parameters,
                                                    experiment_index=index,
                                                    parameter_index=experiment_indices[index].keys()[0])
            try:
                changes["parameter value"] = np.append(changes["parameter value"], parameter_list["parameter value"])
            except KeyError:
                changes["parameter value"] = parameter_list["parameter value"]
            try:
                changes["parameter name"].append(parameter_list["parameter name"])
            except KeyError:
                changes["parameter name"] = [parameter_list["parameter name"]]
            try:
                changes["experiment index"].append(parameter_list["experiment index"])
            except KeyError:
                changes["experiment index"] = [parameter_list["experiment index"]]
            try:
                changes["parameter change"] = np.append(changes["parameter change"], parameter_list["parameter change"])
            except KeyError:
                changes["parameter change"] = parameter_list["parameter change"]
            # changes.append(get_changed_parameters(original_parameters=original,
            #                                       changed_parameters=parameters,
            #                                       experiment_index=index,
            #                                       parameter_index=experiment_indices[index].keys()[0]))
            parameter_id.append(experiment_indices[index].keys()[0])
            parameter_value.append(parameters[index][parameter_id[iter]])
            parameter_change.append(changes["parameter change"][iter])
            ssid_changes.append(list(ssid[index, :]))

        dataset_value_array[data_index[0], :] = np.hstack(data[:])
        # all_parameter_change[data_index[0], :] = np.hstack(changes[:])
        all_parameter_change.append(changes)
        all_parameter_ids.append(parameter_id)
        all_parameter_values.append(parameter_value)
        all_parameter_changes.append(parameter_change)
        all_ssid.append(np.array(ssid_changes))

    dataset_keys = ['values', 'details', 'parameter_ids',
                    'parameter_values', 'parameter_changes', 'ssid',
                    'experiments', 'dataset_id']
    experimental_data = dict(zip(dataset_keys, [dataset_value_array, all_parameter_change, all_parameter_ids,
                                                all_parameter_values, all_parameter_changes, all_ssid,
                                                experiment_index_iter, choose]))
    return experimental_data


def arrange_experimental_data(xss, fss, perturbation_details, experiments_per_set,
                              flux_id=np.array([0, 1, 2, 3, 4, 5]), combination_choice=(), experiment_choice=()):
    """get several data set combinations and
    get data for setting all experimental details for a given combination

    parameters:
    combination_choice - indices of combinations to choose from
    experiment_choice - indices of experiments to choose from to form combinations"""
    number_of_samples = len(xss)
    # get combinations just based on number of experiments in each sample
    data_combinations, data_combination_boolean, \
    experiment_choice, combination_choice = get_data_combinations(xss[0], experiments_per_set,
                                                                  experiment_choice, combination_choice)

    # get values for each combination determined from above permutation
    all_sample_experimental_data = []
    for i_sample in range(0, number_of_samples):
        experimental_data = data_for_each_sample(perturbation_details, experiments_per_set,
                                                 data_combinations, xss[i_sample], fss[i_sample],
                                                 flux_id, combination_choice)
        experimental_data["boolean"] = data_combination_boolean
        all_sample_experimental_data.append(experimental_data)

    return all_sample_experimental_data, experiment_choice, combination_choice

