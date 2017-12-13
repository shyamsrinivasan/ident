import numpy as np
import matplotlib.pyplot as plt
from kotte_model import ident_parameter_name
from kotte_model import flux_based_id


def get_flux_parameter_plot_data(original_data, file_destination=()):
    """calculate and plot the number of data identifying each parameter in each flux"""
    # get parameters identified by each original data set and each combination
    all_boolean_p_id = []
    for len_pos, i_list in enumerate(original_data):
        for i_data in i_list:
            boolean_p_id = [True if j_p in i_data["parameter_ids"] else False for j_p in range(0, 12)]
            all_boolean_p_id.append(boolean_p_id)
    all_boolean_p_id = list(np.transpose(np.array(all_boolean_p_id)))
    # get total data identifying each parameter
    all_boolean_p_id = [sum(j_list) for j_list in all_boolean_p_id]
    # get flux and parameter name for k_p
    flux_name, pid, _ = flux_based_id(range(0, 12))
    parameter_name = ident_parameter_name(pid, flux_name)

    # arrange all data in flux-based dictionary
    unique_flux_names = list(set(flux_name))
    number_of_fluxes = len(unique_flux_names)
    lst_data = [[] for _ in range(number_of_fluxes)]
    for pos, iflux in enumerate(unique_flux_names):
        # get all parameters/fluxes that match iflux in fname
        boolean_id = [j for j, val in enumerate([True if iflux == j_flux else False
                                                 for j_flux in flux_name]) if val]
        # get parameter names for each parameter for each flux in boolean_id (same flux)
        # collect data on the basis of unique fluxes
        lst_data[pos] = {'names':[parameter_name[id] for id in boolean_id],
                         'data':[all_boolean_p_id[id] for id in boolean_id]}

    # plot data for each flux in a separate subplot i.e. number_of_subplots = number_of_fluxes
    nrows = 3
    if number_of_fluxes%3 == 0:
        ncolumns = number_of_fluxes/3
    elif number_of_fluxes%3 == 1:
        ncolumns = (number_of_fluxes - 1)/3
    elif number_of_fluxes%3 == 2:
        ncolumns = (number_of_fluxes + 1) / 3
    else:
        ncolumns = 2
    f, axarr = plt.subplots(nrows, ncolumns, sharex='col')
    total_plots = 0
    for iplot, axis_obj in enumerate(axarr):
        relevant_dict = lst_data[iplot]
        x_data = np.array(relevant_dict["data"])
        y_data = np.arange(len(relevant_dict["names"]))
        axis_obj.barh(y_data, x_data, align='center', color='green', ecolor='black')
        axis_obj.set_yticks(y_data)
        axis_obj.set_yticklabels(relevant_dict["names"])
        axis_obj.invert_yaxis()
        # axis_obj.set_xlabel('Number of Identifying Data Sets')
        axis_obj.set_title(unique_flux_names[iplot])
        total_plots += 1
    axarr[-1].set_xlabel('Number of data sets')
    plt.setp([a.get_xticklabels() for a in axarr[0:-2]], visible=False)
    plt.show()
    if file_destination:
        # save figure to file as png and eps
        plt.savefig(file_destination+'.eps', format='png', dpi=2000)
        plt.savefig(file_destination+'.png', format='eps', dpi=2000)
    return None


def plot_efficient_data():
    return None


def plot_useful_experiments(ident_details, data_list):
    """get most and least useful experiments based on identifiable and non-identifiable datasets"""
    return None


def plot_identifiable_parameter(max_parameter):
    ident_data = max_parameter["data"]
    number_parameters = len(ident_data)
    y_pos = np.arange(number_parameters)
    # x_data = ident_data

    plt.rcdefaults()
    fig, ax = plt.subplots()
    ax.barh(y_pos, ident_data, align='center', color='green', ecolor='black')

    parameter_id = [ident_parameter_name(i) for i in range(0, number_parameters)]

    ax.set_yticks(y_pos)
    ax.set_yticklabels(parameter_id)
    ax.invert_yaxis()
    ax.set_xlabel('Number of Identifying Data Sets')
    ax.set_title('Parameter Identifiability')

    #max_ident_value = max_parameter["maximum"]
    #max_ident_data_id = max_parameter["id"]
    # number of datasets identifying maximum parameters
    #unique_parameter_number_identified = np.unique(np.array())

    #number_max_parameter_data = len(max_ident_data_id)
    return None

