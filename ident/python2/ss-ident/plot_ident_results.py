import numpy as np
import matplotlib.pyplot as plt
from kotte_model import ident_parameter_name
from kotte_model import kotte_experiment_type_name
plt.ion()


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


def plot_on_axis_object_polar(axis_obj, x_data, y_data, data_label, fill_color='b'):
    """plot on polar axis object. does not work on subplots due to the use of plt.methods toplace axis ticks"""
    number_of_experiment_types = len(x_data)
    # get angles for each plot
    angles = [n / float(number_of_experiment_types) * 2 * np.pi for n in range(number_of_experiment_types)]
    # close the circle by adding the first angle at the end
    angles += angles[:1]

    # If you want the first axis to be on top:
    axis_obj.set_theta_offset(np.pi / 2)
    axis_obj.set_theta_direction(-1)

    # set x-ticks (experiment types)
    plt.xticks(angles[:-1], x_data)

    # Draw ylabels
    axis_obj.set_rlabel_position(0)
    max_y_data = max(y_data)
    if max_y_data % 25 != 0:
        if (max_y_data + max_y_data % 25) % 25 != 0:
            max_y_data = (max_y_data - max_y_data % 25) + 25
        else:
            max_y_data = max_y_data + max_y_data % 25

    plt.ylim(0, max_y_data)
    y_ticks = range(25, int(max_y_data)+25, 25)
    y_tick_labels = [str(value) for value in y_ticks]
    plt.yticks(y_ticks, y_tick_labels, color="grey", size=7)

    # set y_data to be similar to angles (close the circle)
    y_data += y_data[:1]

    # plot spider data
    axis_obj.plot(angles, y_data, linewidth=1, linestyle='solid', label=data_label, color=fill_color)
    axis_obj.fill(angles, y_data, fill_color, alpha=0.1)
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
            flux_name = ["flux{}".format(i_flux_info["total"]["flux id"])] * len(y_data)
            flux_choice_id = [i_flux_info["total"]["flux choice"]] * len(y_data)
            parameter_name = ident_parameter_name(y_data, flux_name=flux_name,
                                                  flux_choice_id=flux_choice_id)
            # plot and annotate using plotting function defined above
            plot_on_axis_object(i_axis_obj, x_data, y_data, x_error, x_percent_mean, x_percent_std, noise)
            # set x-axis limits to maximum number of available data sets
            i_axis_obj.set_xlim(0, i_flux_info["percentage"]["data set size"])
            # set custom tick positions
            # i_axis_obj.xaxis.set_ticks(np.arange(0, i_flux_info["percentage"]["data set size"], 40))
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
            flux_choice_id = [i_flux_info["total"]["flux choice"]] * len(y_data)
            parameter_name = ident_parameter_name(y_data, flux_name=flux_name, flux_choice_id=flux_choice_id)
            # plot and annotate using plotting function defined above
            plot_on_axis_object(axarr, x_data, y_data, x_error, x_percent_mean, x_percent_std, noise)
            # set x-axis limits to maximum number of available data sets
            axarr.set_xlim(0, i_flux_info["percentage"]["data set size"])
            # set custom tick positions
            # axarr.xaxis.set_ticks(np.arange(0, i_flux_info["percentage"]["data set size"], 40))
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
                                                  flux_name="flux{}".format(k_parameter_data[0]["total"]["flux id"]),
                                                  flux_choice_id=k_parameter_data[0]["total"]["flux choice"])
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
    all_sample_all_flux_processed_info = flux_based_experiment_info["processed"]
    for j_flux, j_flux_data in enumerate(all_sample_all_flux_processed_info):
        # number_of_subplots = number_of_parameters_in_flux
        number_of_subplots = 1
        number_of_rows = 1
        for k_parameter, k_parameter_data in enumerate(j_flux_data):
            # get parameter name for figure title
            parameter_name = ident_parameter_name(k_parameter,
                                                  flux_name="flux{}".format(k_parameter_data[0]["total"]["flux id"]),
                                                  flux_choice_id=k_parameter_data[0]["total"]["flux choice"])
            # set figure title to parameter name
            figure_title = "flux {}".format(k_parameter_data[0]["total"]["flux id"]) + " " + parameter_name
            # separate plots for each parameter (no subplots)
            f, axarr = plt.subplots(number_of_rows, number_of_subplots, subplot_kw=dict(projection='polar'),
                                    figsize=(6, 4), dpi=100)
            f.text(.5, .975, figure_title, horizontalalignment='center', verticalalignment='top')
            number_of_experiment_positions = len(k_parameter_data)
            all_fill_colors = ['b', 'g', 'y', 'r']
            # collect and plot data from all positions
            for i_position in range(0, number_of_experiment_positions):
                y_data = k_parameter_data[i_position]["percentage"]["mean"]
                x_data = range(0, len(y_data))
                x_data_label = kotte_experiment_type_name(x_data)
                data_label = "experiment {}".format(i_position)
                plot_on_axis_object_polar(axarr, x_data_label, y_data, data_label, all_fill_colors[i_position])
                # add legend
                plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))
            plt.show()
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
            i_axis_obj.set_title('v{} parameters type {}'.format(i_flux_info["total"]["flux id"],
                                                                 i_flux_info["total"]["flux choice"]))
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
            axarr.set_title('v{} parameters type {}'.format(i_flux_info["total"]["flux id"],
                                                            i_flux_info["total"]["flux choice"]))
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


def line_plots_2d(axis_obj, x_data, y_data, x_label=(), y_label=()):
    axis_obj.plot(x_data, y_data)
    if x_label:
        axis_obj.set_xlabel(x_label)
    if y_label:
        axis_obj.set_ylabel(y_label)
    return None


def plot_dynamic_sim_concentrations(dynamic_data, multiple=0):
    """plot dynamic time course profiles of all concentrations present in y vector.
    parameter multiple = 1 when multiple sets of plots are plotted in a subplot"""
    if multiple:
        number_of_plots = len(dynamic_data)
        number_of_rows = 2
        if number_of_plots % number_of_rows != 0:
            number_of_columns = (number_of_plots + 1) / number_of_rows
        else:
            number_of_columns = number_of_plots / number_of_rows
        f, axarr = plt.subplots(number_of_rows, number_of_columns, sharex='col', sharey='row',
                                figsize=(8, 6), dpi=100, facecolor='w',
                                edgecolor='k')
        i_plot = 0
        for i_row in range(0, number_of_rows):
            for i_column in range(0, number_of_columns):
                x_data = dynamic_data[i_plot]["time"]
                y_data = dynamic_data[i_plot]["y"]
                line_plots_2d(axarr[i_row, i_column], x_data, y_data)
                axarr[i_row, i_column].set_title('Data set {}'.format(i_plot + 1))
                if i_row == 0:
                    plt.setp(axarr[i_row, i_column].get_xticklabels(), visible=False)
                if i_column > 0:
                    plt.setp(axarr[i_row, i_column].get_yticklabels(), visible=False)
                if i_row == number_of_rows - 1:
                    axarr[i_row, i_column].set_xlabel('Time (s)')
                if i_column == 0:
                    axarr[i_row, i_column].set_ylabel('Concentrations (a.u.)')
                i_plot += 1
    else:
        f, axis_obj = plt.subplots(1, 1, figsize=(6, 4), dpi=100, facecolor='w', edgecolor='k')
        x_data = dynamic_data["time"]
        y_data = dynamic_data["y"]
        line_plots_2d(axis_obj, x_data, y_data, x_label='Time (s)', y_label='Concentrations (a.u.)')
        axis_obj.set_title('Dynamic Concentrations')
    plt.show()
    return None


def plot_dynamic_sim_fluxes(dynamic_data, multiple=0):
    """plot dynamic time course profiles of all fluxes present in y vector.
        parameter multiple = 1 when multiple sets of plots are plotted in a subplot"""
    if multiple:
        number_of_plots = len(dynamic_data)
        number_of_rows = 2
        if number_of_plots % number_of_rows != 0:
            number_of_columns = (number_of_plots + 1) / number_of_rows
        else:
            number_of_columns = number_of_plots / number_of_rows
        f, axarr = plt.subplots(number_of_rows, number_of_columns, sharex='col', sharey='row',
                                figsize=(8, 6), dpi=100, facecolor='w',
                                edgecolor='k')
        i_plot = 0
        for i_row in range(0, number_of_rows):
            for i_column in range(0, number_of_columns):
                x_data = dynamic_data[i_plot]["time"]
                y_data = dynamic_data[i_plot]["flux"]
                line_plots_2d(axarr[i_row, i_column], x_data, y_data)
                axarr[i_row, i_column].set_title('Data set {}'.format(i_plot + 1))
                if i_row == 0:
                    plt.setp(axarr[i_row, i_column].get_xticklabels(), visible=False)
                if i_column > 0:
                    plt.setp(axarr[i_row, i_column].get_yticklabels(), visible=False)
                if i_row == number_of_rows - 1:
                    axarr[i_row, i_column].set_xlabel('Time (s)')
                if i_column == 0:
                    axarr[i_row, i_column].set_ylabel('Fluxes (a.u.)')
                i_plot += 1
    else:
        f, axis_obj = plt.subplots(1, 1, figsize=(6, 4), dpi=100, facecolor='w', edgecolor='k')
        x_data = dynamic_data["time"]
        y_data = dynamic_data["flux"]
        line_plots_2d(axis_obj, x_data, y_data, x_label='Time (s)', y_label='Fluxes (a.u.)')
        axis_obj.set_title('Dynamic Fluxes')
    plt.show()
    return None


def plot_dynamic_sims(dynamic_data, multiple=0, concentrations=1, fluxes=0):
    if concentrations:
        plot_dynamic_sim_concentrations(dynamic_data, multiple=multiple)

    if fluxes:
        plot_dynamic_sim_fluxes(dynamic_data, multiple=multiple)
    return None


def plot_parameter_values(parameter_values):
    """get box plot of parameter values (determined values vs true values used to create simulations)"""
    number_of_fluxes = len(parameter_values["processed"])
    number_of_plots = number_of_fluxes
    if number_of_fluxes >= 2:
        number_of_rows = 2
    else:
        number_of_rows = 1
    if number_of_plots % number_of_rows != 0:
        number_of_columns = (number_of_plots + 1) / number_of_rows
    else:
        number_of_columns = number_of_plots / number_of_rows
    f, axarr = plt.subplots(number_of_rows, number_of_columns, figsize=(6, 4), dpi=100, facecolor='w', edgecolor='k')
    for i_flux in range(0, number_of_fluxes):
        all_parameter_info = []
        all_parameter_names = []
        all_parameter_true_value = []
        number_parameters = len(parameter_values["processed"][i_flux])
        for i_parameter_info in parameter_values["processed"][i_flux]:
            all_parameter_info.append(i_parameter_info["sample mean"])
            all_parameter_names.append(i_parameter_info["parameter name"])
            all_parameter_true_value.append(i_parameter_info["true value"])
        try:
            bp = axarr[i_flux].boxplot(all_parameter_info)
        except TypeError:
            bp = axarr.boxplot(all_parameter_info)
        for whiskers in bp["whiskers"]:
            whiskers.set(color='k', linewidth=2)
        for flier in bp['fliers']:
            flier.set(marker='o', color='r', alpha=0.5)

        try:
            axarr[i_flux].scatter(range(1, number_parameters+1), all_parameter_true_value,
                                  color='r', marker='o')
            axarr[i_flux].set_xticklabels(all_parameter_names)
            # Remove top axes and right axes ticks
            axarr[i_flux].get_xaxis().tick_bottom()
            axarr[i_flux].get_yaxis().tick_left()
        except TypeError:
            axarr.scatter(range(1, number_parameters + 1), all_parameter_true_value,
                          color='r', marker='o')
            axarr.set_xticklabels(all_parameter_names)
            # Remove top axes and right axes ticks
            axarr.get_xaxis().tick_bottom()
            axarr.get_yaxis().tick_left()
        # axis_obj.set_title()
    plt.show()

    return None
