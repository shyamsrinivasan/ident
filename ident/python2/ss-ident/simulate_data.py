import numpy as np
import itertools as it
import pandas as pd
from collections import defaultdict


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


def data_for_each_sample(exp_df, sample_name, data_combinations, arranged_df=[]):
    """get simulated experimental data for each noisy sample supplied as input argument"""

    for i_combination_number, i_combination in enumerate(data_combinations):
        data_set_name = 'data_set_{}'.format(i_combination_number)
        chosen_series = []
        for i_exp in i_combination:
            get_series = exp_df[exp_df['experiment_id'] == i_exp].loc[sample_name]
            get_series['data_set_id'] = data_set_name
            chosen_series.append(get_series)
        chosen_df = pd.concat(chosen_series, axis=1)
        arranged_df.append(chosen_df.T)

    return arranged_df


# def data_for_each_sample(exp_df, sample_name, data_combinations, all_data_dict, empty_dict):
#     """get simulated experimental data for each noisy sample supplied as input argument"""
#     for i_combination_number, i_combination in enumerate(data_combinations):
#         chosen_df = exp_df.loc(axis=0)[sample_name, i_combination]
#         chosen_df = chosen_df.reset_index(level='experiment_id')
#         # add new column for data set value
#         data_set_name = 'data_set_{}'.format(i_combination_number)
#         chosen_df['data_set_id'] = pd.Series([data_set_name] * len(chosen_df.index.values),
#                                              index=chosen_df.index)
#         chosen_df = chosen_df.reset_index(level='sample_name')
#         temp_dict = chosen_df.to_dict('records')
#         for i_dict in temp_dict:
#             for k, v in it.chain(empty_dict.items(), i_dict.items()):
#                 all_data_dict[k].append(v)
#     return all_data_dict


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
    # all_data = defaultdict(list)
    # empty_dict = {}

    data_df = []
    for i_sample_id, i_sample in enumerate(sample_ids):
        data_df = data_for_each_sample(exp_df, i_sample, data_combinations, data_df)

    df = pd.concat(data_df)

    # rename and reset existing index
    df.index.rename('sample_name', inplace=True)
    df.reset_index(inplace=True)

    # set new index with sample name and data set id
    index_labels = ['sample_name', 'data_set_id']
    df.set_index(index_labels, inplace=True)

    # sort indices
    df.sort_index(level=index_labels, inplace=True)

    return df, index_labels


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
