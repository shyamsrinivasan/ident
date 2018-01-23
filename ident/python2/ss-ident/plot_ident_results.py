import numpy as np
import matplotlib.pyplot as plt
from kotte_model import ident_parameter_name
from kotte_model import kotte_experiment_type_name


def plot_on_axis_object(axis_obj, x_data, y_data, x_error, x_percent_mean, x_percent_std, noise=0):
    """given axis object plot given data on axis object along with all given annotations"""
    # plot bar graphs for x_data vs y_data with x_error error bars
    if noise:
        axis_obj.barh(y_data, x_data, xerr=x_error, align='center', color='blue', ecolor='black')
    else:
        # no error bars
        axis_obj.barh(y_data, x_data, align='center', color='blue', ecolor='black')

    # annotate percentages onto each bar in the graph
    for j_bar in range(0, len(y_data)):
        if noise:
            x_annotation = "{:.2f} + {:.2f}%".format(x_percent_mean[j_bar],
                                                      x_percent_std[j_bar])
        else:
            x_annotation = "{:.2f}%".format(x_percent_mean[j_bar])
        an1 = axis_obj.annotate("", xy=(x_data[j_bar], y_data[j_bar]), xycoords='data',
                                xytext=(x_data[j_bar], y_data[j_bar]), textcoords='data')
        an2 = axis_obj.annotate(x_annotation, xy=(5, .5), xycoords=an1,
                                xytext=(40, 0), textcoords="offset points", size=16, va="center", ha="center")
    # set y axis ticks
    axis_obj.set_yticks(y_data)
    # set y axis tick labels (parameter names)
    return None


def plot_on_axis_object_vertical(axis_obj, x_data, y_data, y_error, y_percent_mean, y_percent_std, noise=0):
    if noise:
        axis_obj.bar(x_data, y_data, yerr=y_error, align='center', color='blue', ecolor='black')
    else:
        # no error bars
        axis_obj.bar(x_data, y_data, align='center', color='blue', ecolor='black')

    # annotate percentages onto each bar in the graph
    for j_bar in range(0, len(x_data)):
        if noise:
            y_annotation = "{:.2f} + {:.2f}%".format(y_percent_mean[j_bar],
                                                      y_percent_std[j_bar])
        else:
            y_annotation = "{:.2f}%".format(y_percent_mean[j_bar])
        an1 = axis_obj.annotate("", xy=(x_data[j_bar], y_data[j_bar]), xycoords='data',
                                xytext=(x_data[j_bar], y_data[j_bar]), textcoords='data')
        an2 = axis_obj.annotate(y_annotation, xy=(5, .5), xycoords=an1,
                                xytext=(40, 0), textcoords="offset points", size=16, va="center", ha="center")
    # set x-ticks
    axis_obj.set_xticks(x_data)
    return None


def parameter_identifibaility_plot(flux_based_parameter_ident, noise=0):
    """plot parameter identifibaility (number of data combinations identifying
    each parameter in each flux. Paramerers for each flux are plotted in separate subplots"""
    all_sample_all_flux_processed_info = flux_based_parameter_ident["processed"]
    number_of_fluxes = len(all_sample_all_flux_processed_info)
    number_of_subplots = number_of_fluxes
    number_of_columns = 1

    # get figure subplots
    f, axarr = plt.subplots(number_of_subplots, number_of_columns, sharex='col',
                            figsize=(8, 6), dpi=100, facecolor='w', edgecolor='k')
    try:
        for i_flux, i_axis_obj in enumerate(axarr):
            i_flux_info = all_sample_all_flux_processed_info[i_flux]
            x_data = i_flux_info["total"]["mean"]
            x_error = i_flux_info["total"]["std"]
            y_data = np.arange(0, len(x_data))
            x_percent_mean = i_flux_info["percentage"]["mean"]
            x_percent_std = i_flux_info["percentage"]["std"]
            # get parameter id/name for y-axis labels
            flux_name = ["flux{}".format(i_flux_info["total"]["flux id"])]*len(y_data)
            parameter_name = ident_parameter_name(y_data, flux_name=flux_name)
            # plot and annotate using plotting function defined above
            plot_on_axis_object(i_axis_obj, x_data, y_data, x_error, x_percent_mean, x_percent_std, noise)
            # set y-axis tick labels
            i_axis_obj.set_yticklabels(parameter_name)
            # set axis title
            i_axis_obj.set_title('v{} parameters'.format(i_flux_info["total"]["flux id"]))
            # invert y-axis
            i_axis_obj.invert_yaxis()
        # set x-axis label
        axarr[-1].set_xlabel('Number of data combinations used for identification')
        # hide x axis tick labels for all but the last subplot sharing x-axes
        plt.setp([a.get_xticklabels() for a in axarr[0:-2]], visible=False)
    except TypeError:
        for i_flux in range(0, number_of_fluxes):
            i_flux_info = all_sample_all_flux_processed_info[i_flux]
            x_data = i_flux_info["total"]["mean"]
            x_error = i_flux_info["total"]["std"]
            y_data = np.arange(0, len(x_data))
            x_percent_mean = i_flux_info["percentage"]["mean"]
            x_percent_std = i_flux_info["percentage"]["std"]
            # get parameter id/name for y-axis labels
            flux_name = ["flux{}".format(i_flux_info["total"]["flux id"])] * len(y_data)
            parameter_name = ident_parameter_name(y_data, flux_name=flux_name)
            # plot and annotate using plotting function defined above
            plot_on_axis_object(axarr, x_data, y_data, x_error, x_percent_mean, x_percent_std, noise)
            # set y-axis tick labels
            axarr.set_yticklabels(parameter_name)
            # set axis title
            axarr.set_title('v{} parameters'.format(i_flux_info["total"]["flux id"]))
            # invert y-axis
            axarr.invert_yaxis()
        # set x-axis label
        axarr.set_xlabel('Number of data combinations used for identification')
    plt.show()
    return None


def parameter_experiment_info_plot(flux_based_experiment_info, noise=0):
    """plot position based contribution from each experiment towards
    identifiable data combinations for each parameter for each flux"""
    all_sample_all_flux_processed_info = flux_based_experiment_info["processed"]
    number_of_fluxes = len(all_sample_all_flux_processed_info)
    for j_flux, j_flux_data in enumerate(all_sample_all_flux_processed_info):
        number_of_parameters_in_flux = len(j_flux_data)
        for k_parameter, k_parameter_data in enumerate(j_flux_data):
            number_of_experiment_positions = len(k_parameter_data)
            number_of_subplots = number_of_experiment_positions
            number_of_rows = 1
            f, axarr = plt.subplots(number_of_rows, number_of_subplots, sharey='row',
                                    figsize=(8, 6), dpi=100, facecolor='w', edgecolor='k')
            # get parameter name for figure title
            parameter_name = ident_parameter_name(k_parameter,
                                                  flux_name="flux{}".format(k_parameter_data[0]["total"]["flux id"]))
            # set figure title to parameter name
            figure_title = "flux {}".format(k_parameter_data[0]["total"]["flux id"]) + " " + parameter_name
            f.text(.5, .975, figure_title, horizontalalignment='center', verticalalignment='top')
            try:
                for i_position, i_axis_obj in enumerate(axarr):
                    x_data = k_parameter_data[i_position]["total"]["mean"]
                    y_data = np.arange(0, len(x_data))
                    x_error = k_parameter_data[i_position]["total"]["std"]
                    x_percent_mean = k_parameter_data[i_position]["percentage"]["mean"]
                    x_percent_error = k_parameter_data[i_position]["percentage"]["std"]
                    # get y-axis labels (experiment types)
                    y_tick_labels = kotte_experiment_type_name(y_data)
                    # plot and annotate using plotting function defined above
                    plot_on_axis_object(i_axis_obj, x_data, y_data, x_error, x_percent_mean, x_percent_error, noise)
                    # set axis title
                    i_axis_obj.set_title('experiment {}'.format(i_position + 1))
                # set x-axis label
                axarr[-1].set_xlabel('Frequency of Experiment Appearance')
                # set y-axis tick label
                axarr[0].set_yticklabels(y_tick_labels)
                # invert y-axis
                axarr[0].invert_yaxis()
            except TypeError:
                for i_position in range(0, number_of_experiment_positions):
                    x_data = k_parameter_data[i_position]["total"]["mean"]
                    y_data = np.arange(0, len(x_data))
                    x_error = k_parameter_data[i_position]["total"]["std"]
                    x_percent_mean = k_parameter_data[i_position]["percentage"]["mean"]
                    x_percent_error = k_parameter_data[i_position]["percentage"]["std"]
                    # plot and annotate using plotting function defined above
                    plot_on_axis_object(axarr, x_data, y_data, x_error, x_percent_mean, x_percent_error, noise)
                    # set axis title
                    axarr.set_title('experiment {}'.format(i_position + 1))
                # set x-axis label
                axarr.set_xlabel('Number of data combinations used for identification')
                # set y-axis tick label
                axarr[0].set_yticklabels(y_tick_labels)
                # invert y-axis
                axarr.invert_yaxis()
    plt.show()
    return None


def parameter_experiment_info_spider(flux_based_experiment_info, noise=0):
    """get spider plots for frequency of different types contributing towards
    identification of different parameters of each flux"""

    return None


def data_utility_plot(data_list, noise=0):
    """collect processed data from all samples for each individual flux and
    plot the utility of the same data combination for all parameters of
    different fluxes using the same data combination"""
    # processed data from all samples
    processed_data = data_list["processed"]
    number_of_fluxes = len(processed_data)
    number_of_subplots = number_of_fluxes
    number_of_columns = 1
    f, axarr = plt.subplots(number_of_subplots, number_of_columns, sharex='col',
                            figsize=(8, 6), dpi=100, facecolor='w', edgecolor='k')
    try:
        max_x_data = []
        # loop through each flux and plot
        for i_flux, i_axis_obj in enumerate(axarr):
            i_flux_info = processed_data[i_flux]
            x_data = i_flux_info["total"]["number"]
            y_data = i_flux_info["total"]["mean"]
            y_error = i_flux_info["total"]["std"]
            y_percent_mean = i_flux_info["percentage"]["mean"]
            y_percent_std = i_flux_info["percentage"]["std"]
            plot_on_axis_object_vertical(i_axis_obj, x_data, y_data, y_error, y_percent_mean, y_percent_std, noise)
            # set axis title
            i_axis_obj.set_title('v{} parameters'.format(i_flux_info["total"]["flux id"]))
            # set x-axis ticks
            max_x_data = max(max_x_data, x_data)
        # set x-axis and y-axis labels
        axarr[-1].set_xlabel('Number of parameters identified')
        axarr[-1].set_ylabel('Number of data combinations')
        # set x-axis ticks
        axarr[-1].set_xticks(max_x_data)
        # set x-axis tick labels
        axarr[-1].set_xticklabels(max_x_data)
    except TypeError:
        max_x_data = []
        for i_flux in range(0, number_of_fluxes):
            i_flux_info = processed_data[i_flux]
            x_data = i_flux_info["total"]["number"]
            y_data = i_flux_info["total"]["mean"]
            y_error = i_flux_info["total"]["std"]
            y_percent_mean = i_flux_info["percentage"]["mean"]
            y_percent_std = i_flux_info["percentage"]["std"]
            plot_on_axis_object_vertical(axarr, x_data, y_data, y_error, y_percent_mean, y_percent_std, noise)
            # set axis title
            axarr.set_title('v{} parameters'.format(i_flux_info["total"]["flux id"]))
            # set x-axis ticks
            max_x_data = max(max_x_data, x_data)
        axarr.set_xlabel('Number of parameters identified')
        axarr.set_ylabel('Number of data combinations')
        # set x-axis ticks
        axarr.set_xticks(max_x_data)
        # set x-axis and y-axis labels
        axarr.set_xticklabels(max_x_data)
    plt.show()
    return None
