import numpy as np
import itertools as it
import pandas as pd
from collections import defaultdict


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


def get_data_combinations(experiment_id, experiments_per_set, experiment_choice, combination_choice):
    """get combinations of data sets based on all perturbations performed on the model"""

    # select only experiment indices provided in experiment_choice
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
    data_combination_boolean = [[True if j_experiment_id in list_1 else False
                                 for j_experiment_id in experiment_id]
                                for list_1 in data_combinations]
    data_combination_boolean = np.array(data_combination_boolean)
    return data_combinations, data_combination_boolean, combination_choice


def data_for_each_sample(exp_df, experiments_per_set, data_combinations, choose):
    """get simulated experimental data for each noisy sample supplied as input argument"""
    new_column_ids = ['experiment_{}_id'.format(i_experiment) for i_experiment in range(0, experiments_per_set)]
    new_exp_df = exp_df.reset_index(level='experiment_id')
    all_column_names = list(new_exp_df.columns.values)
    number_columns = len(all_column_names)
    data_dict = dict(zip(all_column_names, [[] for _ in range(0, number_columns)]))
    data_dict["data_set_id"] = []
    for i_new_column in new_column_ids:
        data_dict[i_new_column] = []
    all_data_sets = []
    for choice_index, data_index in zip(choose, enumerate(data_combinations)):
        data_set_name = 'data_set_{}'.format(data_index[0])
        for i_experiment, i_experiment_index in enumerate(data_index[1]):
            relevant_df = new_exp_df.loc[new_exp_df["experiment_id"] == i_experiment_index, :]
            for i_column in all_column_names:
                data_dict[i_column].append([i_value for i_value in relevant_df[i_column].values][0])
            for j_new_column_id, j_new_column in enumerate(new_column_ids):
                data_dict[j_new_column].append(data_index[1][j_new_column_id])
            data_dict["data_set_id"].append(data_set_name)
        # chosen_experiment_data = exp_df.loc[data_index[1], :]
        # chosen_experiment_data["data_set_id"] = pd.Series([data_set_name] * len(chosen_experiment_data["pep"]),
        #                                                   index=chosen_experiment_data.index)
        # new_data_df = chosen_experiment_data.reset_index(level='experiment_id')
        # experiment_id = list(data_index[1])
        # experiment_id.insert(0, data_set_name)
        # new_index_tuple.append(tuple(experiment_id))
        # new_data_dict = new_data_df.to_dict('list')


    # new_index_tuple = [(i_data, i_exp_1, i_exp_2, i_exp_3) for i_data, i_exp_1, i_exp_2, i_exp_3 in
    #                    zip(data_dict["data_set_id"], data_dict["experiment_0_id"], data_dict["experiment_1_id"],
    #                    data_dict["experiment_2_id"])]
    # index = pd.MultiIndex.from_tuples(new_index_tuple,
    #                                   names=['data_set_id', 'experiment_0_id', 'experiment_1_id', 'experiment_2_id'])
    # del data_dict["data_set_id"]
    # del data_dict["experiment_0_id"]
    # del data_dict["experiment_1_id"]
    # del data_dict["experiment_2_id"]

    # df = pd.DataFrame(data_dict, index=index, columns=data_dict.keys())

    return data_dict


def arrange_experimental_data(exp_df, experiments_per_set, combination_choice=(), experiment_choice=()):
    """get several data set combinations and
    get data for setting all experimental details for a given combination

    parameters:
    combination_choice - indices of combinations to choose from
    experiment_choice - indices of experiments to choose from to form combinations"""
    sample_ids = list(exp_df.index.levels[0])
    # number_samples = len(sample_ids)
    experiment_ids = list(exp_df.index.levels[1])
    # get combinations just based on number of experiments in each sample
    data_combinations, data_combination_boolean, combination_choice = \
        get_data_combinations(experiment_ids, experiments_per_set, experiment_choice, combination_choice)

    # get values for each combination determined from above permutation
    all_data = defaultdict(list)
    d3 = {}

    for i_sample_id, i_sample in enumerate(sample_ids):
        for i_combination_number, i_combination in enumerate(data_combinations):
            chosen_df = exp_df.loc(axis=0)[i_sample, i_combination]
            chosen_df = chosen_df.reset_index(level='experiment_id')
            # add new column for data set value
            data_set_name = 'data_set_{}'.format(i_combination_number)
            chosen_df['data_set_id'] = pd.Series([data_set_name] * len(chosen_df.index.values),
                                                             index=chosen_df.index)
            chosen_df = chosen_df.reset_index(level='sample_name')
            temp_dict = chosen_df.to_dict('records')
            for i_dict in temp_dict:
                for k, v in it.chain(d3.items(), i_dict.items()):
                    all_data[k].append(v)
    # create multi index tuples
    df_tuples = [(j_sample, j_data_set, j_experiment) for j_sample, j_data_set, j_experiment in
                 zip(all_data["sample_name"], all_data["data_set_id"], all_data["experiment_id"])]
    index_labels = ['sample_name', 'data_set_id', 'experiment_id']
    # create multi index
    index = pd.MultiIndex.from_tuples(df_tuples, names=index_labels)
    # remove redundant columns
    del all_data["sample_name"]
    del all_data["data_set_id"]
    del all_data["experiment_id"]
    # create data frame
    df = pd.DataFrame(all_data, index=index, columns=all_data.keys())

    return all_sample_experimental_data, experiment_choice, combination_choice


def data_for_each_sample_numerical(perturbation_details, experiments_per_set,
                         data_combinations, xss, fss, flux_id, choose):
    """return data required for numerical determination of identifiability"""
    # collect data for combinations in data_combinations for each sample in xss
    parameters = perturbation_details["values"]

    all_data_values = []
    experiment_index_iter = []
    for choice_index, data_index in zip(choose, enumerate(data_combinations)):
        experiment_index_iter.append(data_index[1])
        data = []
        for iter, index in enumerate(data_index[1]):
            data.append(np.hstack((parameters[index]["ac"],
                                   xss[index],
                                   fss[index][flux_id])))
        all_data_values.append(data)
    dataset_keys = ['values', 'dataset_id', 'experiments']
    experimental_data = dict(zip(dataset_keys, [all_data_values, choose, experiment_index_iter]))
    return experimental_data


def arrange_experimental_data_numerical(xss, fss, perturbation_details, experiments_per_set,
                              flux_id=np.array([0, 1, 2, 3, 4, 5]), combination_choice=(), experiment_choice=()):
    """arrange experimental is a list with len(list) = number of experiments needed for identifying each flux"""
    number_of_samples = len(xss)
    # get combinations just based on number of experiments in each sample
    data_combinations, data_combination_boolean, \
    experiment_choice, combination_choice = get_data_combinations(xss[0], experiments_per_set,
                                                                  experiment_choice, combination_choice)

    # get values for each combination determined from above permutation
    all_sample_experimental_data = []
    for i_sample in range(0, number_of_samples):
        experimental_data = data_for_each_sample_numerical(perturbation_details, experiments_per_set,
                                                           data_combinations, xss[i_sample], fss[i_sample],
                                                           flux_id, combination_choice)
        experimental_data["boolean"] = data_combination_boolean
        all_sample_experimental_data.append(experimental_data)
    return all_sample_experimental_data, experiment_choice, combination_choice

