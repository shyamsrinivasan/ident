from create_experiment_data import retrieve_experimental_data_from_file
from collections import defaultdict
import pandas as pd
import itertools as it


def get_original_experiments(exp_file_name):
    # original set of experiments
    exp_df = retrieve_experimental_data_from_file(data_file_name=exp_file_name,
                                                  multi_index_label=['sample_name', 'experiment_id'])
    all_possible_experiments = exp_df.index.levels[1].tolist()
    return all_possible_experiments


def get_ident_experiments(ident_df, all_experiment_ids, max_number_experiments, all_parameter_info, all_parameter_flux,
                          parameter_names):
    """get occurrence of each experiment (not experiment type) in identifying data
    for each parameter in a given data set"""
    all_parameter_names = ident_df["parameter_name"].unique().tolist()

    # set up list of experiment positions and experiment ids
    all_exp_column_ids = ['experiment_{}_id'.format(i_experiment) for i_experiment in range(0, max_number_experiments)]

    # collect experiment occurrence info for creating data frames
    empty_dict = {}
    # all_parameter_info = []
    # all_parameter_flux = []
    for i_parameter, i_parameter_name in enumerate(all_parameter_names):
        # get all data sets identifying each parameter
        all_ident_exp_info = defaultdict(list)
        identifying_df = ident_df[(ident_df["parameter_name"] == i_parameter_name) & (ident_df["identified"])]
        for i_expeirment, i_experiment_pos_name in enumerate(all_exp_column_ids):
            # get frequency of each experiment at i_experiment position
            try:
                exp_frequency = identifying_df[i_experiment_pos_name].value_counts()
            except KeyError:
                # if experiment_pos_name is not available in df
                # (only 2 experiments against maximum of 3 are required)
                # create a series with all zero values
                exp_frequency = pd.Series([0] * len(all_experiment_ids), index=all_experiment_ids,
                                          name=i_experiment_pos_name)
            # add non-existent experiments to series with zero values
            available_experiments = exp_frequency.index.tolist()
            missing_experiments = [i_exp for i_exp in all_experiment_ids if i_exp not in available_experiments]
            missing_experiment_value = [0] * len(missing_experiments)
            missing_series = pd.Series(missing_experiment_value, index=missing_experiments)
            new_experiment_frequency = exp_frequency.append(missing_series).rename(exp_frequency.name)
            # sort series index before writing into dictionary
            new_experiment_frequency.sort_index(inplace=True)
            experiment_ids = new_experiment_frequency.index.tolist()
            experiment_frequency = [int(i_value) for i_value in new_experiment_frequency.values.tolist()]
            experiment_pos = [i_experiment_pos_name] * len(experiment_ids)
            new_dict = dict(zip(['experiment_id', 'frequency', 'experiment_pos_id'],
                                [experiment_ids, experiment_frequency, experiment_pos]))
            for key, value in it.chain(empty_dict.items(), new_dict.items()):
                for i_value in value:
                    all_ident_exp_info[key].append(i_value)
        all_parameter_info.append(all_ident_exp_info)
        all_parameter_flux.append(identifying_df["flux_name"].unique()[0])
        parameter_names.append(i_parameter_name)
    return all_parameter_info, all_parameter_flux, parameter_names


def exp_design_info(list_of_files, original_experiment_file, write_to_file_name, max_number_experiments):
    """get df of all experiments contributing to identification of all parameters in all fluxes
    (ident info for each flux is given in a separate file)"""

    # get all info on original experiments
    all_experiment_ids = get_original_experiments(original_experiment_file)

    # get ident details
    ident_index_label = ['sample_name', 'data_set_id']
    all_parameter_info = []
    all_parameter_flux = []
    all_parameter_names = []
    for i_file in list_of_files:
        # retrieve ident info from file
        ident_df = retrieve_experimental_data_from_file(i_file, ident_index_label)

        # get experiment contribution for each parameter in each flux
        all_parameter_info, all_parameter_flux, all_parameter_names = \
            get_ident_experiments(ident_df, all_experiment_ids, max_number_experiments=max_number_experiments,
                                  all_parameter_info=all_parameter_info, all_parameter_flux=all_parameter_flux,
                                  parameter_names=all_parameter_names)
    all_frequency = [i_parameter["frequency"] for i_parameter in all_parameter_info]
    col_ind_tuple = zip(all_parameter_info[0]['experiment_pos_id'], all_parameter_info[0]['experiment_id'])
    col_ind = pd.MultiIndex.from_tuples(col_ind_tuple, names=['experiment_pos_id', 'experiment_id'])
    flux_tuple = zip(all_parameter_flux, all_parameter_names)
    row_ind = pd.MultiIndex.from_tuples(flux_tuple, names=['flux_name', 'parameter_name'])
    df = pd.DataFrame(all_frequency, index=row_ind, columns=col_ind)
    df.to_csv(write_to_file_name, index_label=['flux_name', 'parameter_name'])
    return df
