import numpy as np


def truncate_values(f, n=3):
    """truncates floats to n specified values after the decimal"""
    if not np.isnan(f):
        s = '{}'.format(f)  # convert float to string
        if 'e' in s or 'E' in s:
            return float('{0:.{1}f}'.format(f, n))
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


def boolean_ident_info(ident_values, number_of_parameters):
    """get absolute values for each identifiability function and convert them to boolean arrays
    to denote whether a given combination of experimental data can identify a given parameter"""
    signed_ident_values = np.sign(ident_values)
    p_list = [[p_id for p_id, val in enumerate(data_set) if val > 0] for data_set in signed_ident_values]
    p_list_boolean = [[True if parameter_id in list_1 else False for parameter_id in range(0, number_of_parameters)]
                      for list_1 in p_list]
    return p_list, np.array(p_list_boolean)


def run_flux_ident(ident_function_list, data, flux_id=(), flux_choice=()):
    """test identifibaility using each function in function list for data set in set"""
    number_of_parameters = 0
    try:
        number_of_parameters_per_flux = [None] * len(ident_function_list)
    except TypeError:
        number_of_parameters_per_flux = [None]
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
            number_of_expressions = len(ident_value[0])
            number_of_parameters_per_flux[iterator] = len(ident_value)
            iterator += 1
            number_of_parameters += len(ident_value)
    except TypeError:
        ident_value = ident_function_list(data)
        ident_value_list.append(ident_value)
        flux_id_list.append(flux_id)
        flux_choice_list.append(flux_choice)
        number_of_expressions = len(ident_value[0])
        number_of_parameters_per_flux[iterator] = len(ident_value)
        number_of_parameters += len(ident_value)

    ident_value_array = np.zeros((number_of_parameters, number_of_expressions))
    irow = 0
    all_flux_ident = []
    for iflux in ident_value_list:
        truncated_ident_value = call_truncate_method(iflux, len(iflux))
        all_flux_ident.append(truncated_ident_value)
        nrows, ncolumns = np.shape(truncated_ident_value)
        ident_value_array[irow:(irow + nrows), :] = truncated_ident_value
        irow += nrows
    return ident_value_array, number_of_parameters_per_flux, ncolumns, all_flux_ident, \
           flux_id_list, flux_choice_list


def get_ident_value(ident_function_list, experimental_data_list, original_data_set_id, flux_ids, flux_choice):
    """perform idetifibaility for each data set by looping through data sets and
    using them to test each function in identifibaility list.
    Also process results into numpy arrays for future parsing, processing and analysis"""
    all_data_ident_lists = []
    all_data_all_fun_ident_value = []
    try:
        number_of_parameters_per_flux = [0]*len(ident_function_list)
    except TypeError:
        number_of_parameters_per_flux = [0]
    number_of_expressions_per_parameter = 0
    for index, data_set in enumerate(experimental_data_list):
        print('Identifiability for Dataset {} of {}: Original ID: {}\n'.format(index + 1,
                                                                               len(experimental_data_list),
                                                                               original_data_set_id[index]))
        identifiability_values, number_of_parameters_per_flux, \
        number_of_expressions_per_parameter, all_fun_ident_values, \
        flux_id_list, flux_choice_list = run_flux_ident(ident_function_list, data_set, flux_ids, flux_choice)
        all_data_ident_lists.append(identifiability_values)
        all_data_all_fun_ident_value.append(all_fun_ident_values)

    # classify ident data based on number of ident functions instead of by data sets
    # convert to multi dimensional list = data set size x number of parameters per flux x number of fluxes/functions
    number_of_data_sets = len(experimental_data_list)
    # all_fun_ident_list = []
    all_fun_array_list = []
    try:
        for ifun in range(0, len(ident_function_list)):
            all_ident_data_in_fun = []
            for idata in all_data_all_fun_ident_value:
                all_ident_data_in_fun.append(idata[ifun])
            # all_fun_ident_list.append(all_ident_data_in_fun)
            # array - number of data sets x number of parameters per identifiability function/flux
            ident_fun_array = np.zeros((number_of_data_sets, number_of_parameters_per_flux[ifun], 3))
            for idata_id, idata_set_value in enumerate(all_ident_data_in_fun):
                ident_fun_array[idata_id, :, 0] = idata_set_value[:, 0]
                ident_fun_array[idata_id, :, 1] = idata_set_value[:, 1]
                ident_fun_array[idata_id, :, 2] = idata_set_value[:, 2]
            all_fun_array_list.append(ident_fun_array)
    except TypeError:
        ifun = 0
        all_ident_data_in_fun = []
        for idata in all_data_all_fun_ident_value:
            all_ident_data_in_fun.append(idata[ifun])
        # all_fun_ident_list.append(all_ident_data_in_fun)
        # array - number of data sets x number of parameters per identifiability function/flux
        ident_fun_array = np.zeros((number_of_data_sets, number_of_parameters_per_flux[ifun], 3))
        for idata_id, idata_set_value in enumerate(all_ident_data_in_fun):
            ident_fun_array[idata_id, :, 0] = idata_set_value[:, 0]
            ident_fun_array[idata_id, :, 1] = idata_set_value[:, 1]
            ident_fun_array[idata_id, :, 2] = idata_set_value[:, 2]
        all_fun_array_list.append(ident_fun_array)


    total_parameters_identified = sum(number_of_parameters_per_flux)
    all_identifiability_values = \
        np.zeros((number_of_data_sets * total_parameters_identified, number_of_expressions_per_parameter))
    array_index = 0
    for idata in all_data_ident_lists:
        all_identifiability_values[array_index:(array_index + total_parameters_identified), :] = idata
        array_index += total_parameters_identified
    return all_identifiability_values, number_of_parameters_per_flux, all_fun_array_list


def one_sample_ident_fun(ident_fun_list, sample_data, choose, flux_ids, flux_choice):
    """get identifibaility information for a single sample of experimental data
    using identifiablity function list passed as input args"""
    if choose:
        try:
            chosen_values = list(sample_data["values"][:choose, :])
        except TypeError:
            iter_chosen_value = []
            for indexes in range(0, len(choose)):
                iter_chosen_value.append(list(sample_data["values"][indexes, :]))
            chosen_values = iter_chosen_value[:]
    else:
        chosen_values = list(sample_data["values"][:, :])
    # run identification function through every chosen data combination supplied as input
    _, parameters_per_flux, all_fun_array_list = get_ident_value(ident_fun_list, chosen_values, choose,
                                                                 flux_ids, flux_choice)
    return parameters_per_flux, all_fun_array_list


def multi_sample_ident_fun(ident_fun_list, all_data, choose, flux_ids, flux_choice):
    """perform identifibaility analysis for multiple samples by
    looping through each experimental data sample for a list of identifibaility functions"""
    all_sample_ident_details = []
    for i_sample, sample_data in enumerate(all_data):
        parameters_per_flux, one_sample_all_fun_array_list = one_sample_ident_fun(ident_fun_list, sample_data,
                                                                                  choose, flux_ids, flux_choice)
        # process info from each sample (number of samples > 1 when noisy data is used) of experimental data sets
        all_fun_ident_details = []
        for ifun, ifun_data in enumerate(one_sample_all_fun_array_list):
            number_of_data_sets, number_of_parameters, _ = ifun_data.shape
            # process denominator info only
            _, plist_boolean = boolean_ident_info(ifun_data[:, :, -1], number_of_parameters)
            ident_details = {"boolean": plist_boolean,
                             "values": ifun_data,
                             "parameters": number_of_parameters,
                             "flux id": flux_ids[ifun],
                             "flux choice": flux_choice[ifun]}
            all_fun_ident_details.append(ident_details)
        all_sample_ident_details.append(all_fun_ident_details)
    return all_sample_ident_details
