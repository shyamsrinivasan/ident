import numpy as np
import pandas as pd
from collections import defaultdict
import itertools as it
from names_strings import ident_parameter_name


def truncate_values(f, n=3):
    """truncates floats to n specified values after the decimal"""
    if not np.isnan(f):
        if type(f) is not np.ndarray:
            s = '{}'.format(f)  # convert float to string
            if 'e' in s or 'E' in s:
                return float('{0:.{1}f}'.format(f, n))
        else:
            s = '{}'.format(f[0])  # convert np.ndarray to string
            if 'e' in s or 'E' in s:
                return float('{0:.{1}f}'.format(f[0], n))
        i, p, d = s.partition('.')
        return float('.'.join([i, (d+'0'*n)[:n]]))
    else:
        return f


def call_truncate_method(ident_value_list, parameter_count, expression_count=3):
    """calculate truncate_values using map on list of values"""
    flux_ident_value = np.zeros((parameter_count, expression_count))
    for i, j in enumerate(ident_value_list):
        trunc_value = map(truncate_values, j)
        # trunc_value = map(float, trunc_value)
        flux_ident_value[i, :] = np.array(trunc_value)
    return flux_ident_value


def run_flux_ident(ident_function_list, data, flux_id=(), flux_choice=()):
    """test identifibaility using each function in function list for data set in set"""
    ident_value_list = []
    flux_id_list = []
    flux_choice_list = []
    iterator = 0
    if not flux_id:
        flux_id = range(1, len(ident_function_list)+1)
        flux_choice = [0] * len(ident_function_list)
    try:
        for func, i_d in zip(ident_function_list, flux_id):
            ident_value = func(data)
            ident_value_list.append(ident_value)
            flux_id_list.append(i_d)
            flux_choice_list.append(flux_choice[iterator])
            iterator += 1
    except TypeError:
        ident_value = ident_function_list(data)
        ident_value_list.append(ident_value)
        flux_id_list.append(flux_id)
        flux_choice_list.append(flux_choice)

    all_flux_ident = []
    for iflux in ident_value_list:
        truncated_ident_value = call_truncate_method(iflux, len(iflux))
        ident_value_list = [np.array(i_parameter) for i_parameter in list(truncated_ident_value)]
        all_flux_ident.append(ident_value_list)
    return all_flux_ident, flux_id_list, flux_choice_list


def get_ident_value(ident_function_list, experimental_data_list, flux_ids, flux_choice):
    """perform idetifibaility for each data set by looping through data sets and
    using them to test each function in identifibaility list.
    Also process results into numpy arrays for future parsing, processing and analysis"""
    # all_data_sets = len(experimental_data_list)
    number_data_sets = (len(experimental_data_list))
    all_data_ident_lists = []
    # number_of_expressions_per_parameter = 0
    for j_data_set_id, j_data_set in enumerate(experimental_data_list):
        print('Identifiability for Dataset {} of {}: \n'.format(j_data_set_id, number_data_sets))
        identifiability_values, flux_id_list, flux_choice_list = \
            run_flux_ident(ident_function_list, j_data_set, flux_ids, flux_choice)
        all_data_ident_lists.append(identifiability_values)
    number_of_parameters_per_flux = [len(i_flux_info) for i_flux_info in all_data_ident_lists[0]]

    return all_data_ident_lists, number_of_parameters_per_flux


def collect_data(exp_df, j_sample, numerical=0):
    """collect data from df from all data sets in single sample to
    peform identiiability analysis with cas and numerical method"""
    idx = pd.IndexSlice
    all_data_set_ids = exp_df.index.levels[1].tolist()
    # number_data_sets = (len(all_data_sets))
    all_exp_data = []
    if numerical:
        # avoid converting list of lists to single list
        for j_data_set, data_set_id in enumerate(all_data_set_ids):
            ident_data = exp_df.loc[idx[j_sample, data_set_id],
                                    ['acetate', 'pep', 'fdp', 'E', 'v1', 'v2', 'v3', 'v5']].values.tolist()
            single_list = [i_exp_data for i_exp_data in ident_data]
            all_exp_data.append(single_list)
    else:
        for j_data_set, data_set_id in enumerate(all_data_set_ids):
            ident_data = exp_df.loc[idx[j_sample, data_set_id],
                                    ['acetate', 'pep', 'fdp', 'E', 'v1', 'v2', 'v3', 'v5']].values.tolist()
            single_list = [i_variable for i_exp_data in ident_data for i_variable in i_exp_data]
            all_exp_data.append(single_list)
    return all_exp_data, all_data_set_ids


def collect_ident_data(j_sample_name, j_sample_ident_data, flux_ids, flux_choice, all_data_dict, empty_dict):
    """collect identifiability information from processing info from
    using experimental data from all samples"""

    temp_dict = {}
    for j_data_set, j_data_set_info in enumerate(j_sample_ident_data):
        data_set_name = 'data_set_{}'.format(j_data_set)
        for j_flux, j_flux_data in enumerate(j_data_set_info):
            flux_name = 'flux{}'.format(flux_ids[j_flux])
            # all_parameter names
            all_parameter_names = [ident_parameter_name(j_parameter,
                                                        flux_name,
                                                        flux_choice[j_flux])
                                   for j_parameter in range(0, len(j_flux_data))]
            for i_parameter, i_parameter_info in enumerate(j_flux_data):
                temp_dict["flux_name"] = flux_name
                temp_dict["flux_choice"] = flux_choice[j_flux]
                # replace with call to parameter name file
                temp_dict["parameter_name"] = all_parameter_names[i_parameter]
                i_parameter_nr, i_parameter_dr, i_parameter_value = i_parameter_info
                temp_dict["parameter_nr"] = i_parameter_nr
                temp_dict["parameter_dr"] = i_parameter_dr
                temp_dict["parameter_value"] = i_parameter_value
                temp_dict["data_set_id"] = data_set_name
                temp_dict["sample_name"] = j_sample_name
                if i_parameter_value > 0:
                    temp_dict["identified"] = True
                else:
                    temp_dict["identified"] = False
                for key, value in it.chain(empty_dict.items(), temp_dict.items()):
                    all_data_dict[key].append(value)

    return all_data_dict


def multi_sample_ident_fun(ident_fun_list, all_data_df, flux_ids, flux_choice):
    """perform identifibaility analysis for multiple samples by
    looping through each experimental data sample for a list of identifibaility functions"""
    reset_df = all_data_df.reset_index('experiment_id')
    # lexographic ordering of df indices
    reset_df.sort_index(level='data_set_id', inplace=True)
    reset_df.sort_index(level='sample_name', inplace=True)
    sample_ids = list(reset_df.index.levels[0])
    number_samples = (len(sample_ids))
    all_sample_ident_details = []
    for i_sample, i_sample_id in enumerate(sample_ids):
        print('Identifiability analysis with Data Sample Number {} of {}\n'.format(i_sample, number_samples))
        # collect experimental data from all data sets
        all_exp_data, _ = collect_data(reset_df, i_sample_id)
        # run identifiability with i_sample_data
        all_ident_values, _ = get_ident_value(ident_fun_list, all_exp_data, flux_ids, flux_choice)
        all_sample_ident_details.append(all_ident_values)

    # initial information processing to get dictionary of relevant info for each flux and each parameter
    all_data = defaultdict(list)
    empty_dict = {}
    for i_sample, sample_data in enumerate(all_sample_ident_details):
        sample_name = 'sample_{}'.format(i_sample)
        # number_data_sets = len(sample_data)
        all_data = collect_ident_data(sample_name, sample_data, flux_ids, flux_choice, all_data, empty_dict)

    return all_data


def data_numerical_ident(exp_df, sample_name):
    """create data sets of appropriate form for numerical identifiability"""
    exp_df.sort_index(level='data_set_id', inplace=True)
    exp_df.sort_index(level='sample_name', inplace=True)
    all_exp_df, data_set_ids = collect_data(exp_df, j_sample=sample_name, numerical=1)
    all_data_info = {"value": all_exp_df, "data_set_id": data_set_ids}
    return all_data_info
