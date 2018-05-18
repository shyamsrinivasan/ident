import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.font_manager as fnt
from names_strings import ident_parameter_name
from names_strings import kotte_experiment_type_name
from names_strings import variable_name
plt.ion()


def plot_on_axis_object(axis_obj, x_data, y_data, x_error, noise=0):
    """given axis object plot given data as horizontal bar on axis object"""
    # plot bar graphs for x_data vs y_data with x_error error bars
    if noise:
        axis_obj.barh(y_data, x_data, xerr=x_error, align='center', ecolor='black', capsize=0.3,
                      **{"alpha": 0.4, "capstyle": 'projecting'})
    else:
        # no error bars
        axis_obj.barh(y_data, x_data, align='center', height=0.2, color='blue', **{"alpha": 0.4})
    return None


def set_hbar_axis_properties(axis_obj, y_data, y_tick_label, x_max, x_percent_mean, x_label=[], x_percent_std=[], figure_title=[]):
    """set horizontal bar plot properties along with annotations for given """
    # set x-axis limits
    axis_obj.set_xlim(0, x_max)
    # recompute data limits
    axis_obj.relim()
    axis_obj.autoscale_view()
    # set custom tick positions
    # axis_obj.xaxis.set_ticks(np.arange(0, x_max, 40))

    # annotate percentages onto each bar in the graph
    f_p = fnt.FontProperties(size=14, weight='demibold')
    for j_bar, p in enumerate(axis_obj.patches):
        # annotate percentages onto bar graphs
        if x_percent_std:
            y_annotation = "{:.2f} + {:.2f}%".format(x_percent_mean[j_bar],
                                                     x_percent_std[j_bar])
        else:
            y_annotation = "{:.2f}%".format(x_percent_mean[j_bar])
        axis_obj.annotate(y_annotation, (p.get_width() / 2, p.get_y() + p.get_height() / 2),
                          ha='center', va='center', xycoords='data', **{"color": 'black', "font_properties": f_p})
    # set y axis ticks
    axis_obj.set_yticks(y_data)
    # set y axis tick labels (parameter names)
    # y_tick_label_p = fnt.FontProperties(size=14, weight='demibold')
    axis_obj.set_yticklabels(y_tick_label, **{"font_properties": f_p})
    # set x-axis label
    if x_label:
        axis_obj.set_xlabel(x_label[0], **{"font_properties": f_p})
    axis_obj.tick_params(axis='both', labelsize=14, color='grey')
    # set axis title
    if figure_title:
        axis_obj.set_title(figure_title, **{"font_properties": f_p})
    # invert y-axis
    axis_obj.invert_yaxis()

    return None


def plot_on_axis_object_vertical(axis_obj, x_data, y_data, y_error, y_percent_mean, y_percent_std, y_max, noise=0):
    if noise:
        axis_obj.bar(x_data, y_data, width=.5, yerr=y_error, align='center')
    else:
        # no error bars
        axis_obj.bar(x_data, y_data, width=.5, align='center')
    # set y-axis limits
    axis_obj.set_ylim(0, y_max)
    # recompute data limits
    axis_obj.relim()
    axis_obj.autoscale_view()
    # set custom tick positions
    # axis_obj.yaxis.set_ticks(np.arange(0, y_max, 40))

    # annotate percentages onto each bar in the graph
    for j_bar, p in enumerate(axis_obj.patches):
        p.set_facecolor('blue')
        p.set_edgecolor('k')
        # set hatch pattern
        # p.set_hatch('/')

        # annotate percentages onto bar graphs
        if noise:
            y_annotation = "{:.2f} + {:.2f}%".format(y_percent_mean[j_bar],
                                                      y_percent_std[j_bar])
        else:
            y_annotation = "{:.2f}%".format(y_percent_mean[j_bar])
        axis_obj.annotate(y_annotation, (p.get_x() + p.get_width()/2, p.get_height()*1.05),
                          ha='center', va='center', textcoords='data')
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
    axis_obj.set_thetagrids(np.degrees(angles[:-1]), x_data)

    # Draw ylabels
    axis_obj.set_rlabel_position(0)
    y_data_bool = [True if iy <= 0.0 else False for iy in y_data]
    # exception here when all y_data = 0.0

    # set y_data to be similar to angles (close the circle)
    y_data += y_data[:1]

    # plot spider data
    # exception here when all y_data = 0.0
    if not all(y_data_bool):
        axis_obj.plot(angles, y_data, linewidth=1, linestyle='solid', label=data_label, color=fill_color)
        axis_obj.fill(angles, y_data, fill_color, alpha=0.1)
        max_y_data = max(y_data)
        if max_y_data % 10 != 0:
            if (max_y_data + max_y_data % 10) % 10 != 0:
                max_y_data = (max_y_data - max_y_data % 10) + 10
            else:
                max_y_data = max_y_data + max_y_data % 10
    else:
        max_y_data = 0.0
    return max_y_data


def set_polar_axis_limits(axis_object, y_limit):
    """set axis limits for polar plots"""
    axis_object.set_rmax(y_limit)
    y_ticks = range(0, int(y_limit) + 10, 10)
    axis_object.set_rticks(y_ticks)
    axis_object.set_rlabel_position(22.5)
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
            x_max = i_flux_info["percentage"]["data set size"]
            # get parameter id/name for y-axis labels
            flux_name = ["flux{}".format(i_flux_info["total"]["flux id"])] * len(y_data)
            flux_choice_id = [i_flux_info["total"]["flux choice"]] * len(y_data)
            parameter_name = ident_parameter_name(y_data, flux_name=flux_name,
                                                  flux_choice_id=flux_choice_id)
            # plot and annotate using plotting function defined above
            plot_on_axis_object(i_axis_obj, x_data, y_data, x_error, x_percent_mean, x_percent_std, x_max, noise)
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
            x_max = i_flux_info["percentage"]["data set size"]
            # get parameter id/name for y-axis labels
            flux_name = ["flux{}".format(i_flux_info["total"]["flux id"])] * len(y_data)
            flux_choice_id = [i_flux_info["total"]["flux choice"]] * len(y_data)
            parameter_name = ident_parameter_name(y_data, flux_name=flux_name, flux_choice_id=flux_choice_id)
            # plot and annotate using plotting function defined above
            plot_on_axis_object(axarr, x_data, y_data, x_error, x_percent_mean, x_percent_std, x_max, noise)
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


def exp_info_plot(info_dict):
    """plot experiment contribution frequency for each parameter in a polar plot"""
    number_parameters = len(info_dict["names"])
    if number_parameters >= 3:
        number_of_columns = 3
    else:
        number_of_columns = 1
    if number_parameters % number_of_columns != 0:
        number_of_rows = (number_parameters + 1) / number_of_columns
    else:
        number_of_rows = number_parameters / number_of_columns
    # figure = plt.figure(figsize=(6, 4))
    # inner_grid = gridspec.GridSpec(number_of_rows, number_of_columns, wspace=0.2, hspace=0.2)
    f, ax = plt.subplots(number_of_rows, number_of_columns, subplot_kw=dict(projection='polar'),
                         figsize=(6, 4), dpi=100, gridspec_kw={"wspace": 0.2, "hspace": 0.2})
    all_max_y_data = []
    for i_parameter, i_parameter_info in enumerate(info_dict["exp_info"]):
        # ax = plt.Subplot(figure, inner_grid[i_parameter])
        fill_colors = ['b', 'g', 'y', 'r']
        # plot data from all positions
        pos_labels = []
        for i_pos, (i_pos_key, i_pos_val) in enumerate(i_parameter_info.items()):
            pos_labels.append(i_pos_key)
            y_data = i_pos_val["frequency"]
            x_labels = i_pos_val["names"]
            max_y_data = plot_on_axis_object_polar(ax[i_parameter], x_data=x_labels, y_data=y_data,
                                                   data_label=pos_labels[-1], fill_color=fill_colors[i_pos])
            all_max_y_data.append(max_y_data)

    # set axis limits for all polar plots on subplot
    for i_parameter, _ in enumerate(info_dict["exp_info"]):
        set_polar_axis_limits(ax[i_parameter], max(all_max_y_data))

    # add legend
    plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))
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
            y_max = i_flux_info["percentage"]["data set size"]
            plot_on_axis_object_vertical(i_axis_obj, x_data, y_data, y_error,
                                         y_percent_mean, y_percent_std, y_max, noise)
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
            y_max = i_flux_info["percentage"]["data set size"]
            plot_on_axis_object_vertical(axarr, x_data, y_data, y_error,
                                         y_percent_mean, y_percent_std, y_max, noise)
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


def plot_dynamic_sim_concentrations(dynamic_data, multiple=0, noise=0):
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
        if noise:
            # plt only one of many noisy slutions (Figure out a way to plot all)
            line_plots_2d(axis_obj, x_data, y_data[0], x_label='Time (s)', y_label='Concentrations (a.u.)')
        else:
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


def multi_box_plot(figure, axis_object, plot_data, noise=1, plot_positions=[]):
    """generate box plot of distribution of data in plot_data"""
    if not noise:
        axis_object = plt.Subplot(figure, axis_object)
    if plot_positions:
        bp = axis_object.boxplot(plot_data, positions=plot_positions, whis='range')
    else:
        bp = axis_object.boxplot(plot_data)
    for whiskers in bp["whiskers"]:
        whiskers.set(color='k', linewidth=2)
    for flier in bp['fliers']:
        flier.set(marker='o', color='r', alpha=0.5)
    # Remove top axes and right axes ticks
    axis_object.get_xaxis().tick_bottom()
    axis_object.get_yaxis().tick_left()
    figure.add_subplot(axis_object)
    plt.show()
    return axis_object


def get_true_boolean_data(raw_ident_data, raw_ident_boolean_data, data_set_indices_to_plot):
    """get identifying samples for each data set and collate in list of size = number_data_sets"""
    # number_samples, number_data_sets = raw_ident_data.shape
    # go through each data set and pick only values from samples that have true values
    true_ident_data = []  # range(0, number_data_sets)
    selected_data_id = []
    for i_data in data_set_indices_to_plot:
        if i_data < raw_ident_boolean_data.shape[1]:
            true_ident_data.append(raw_ident_data[raw_ident_boolean_data[:, i_data], i_data])
            selected_data_id.append(i_data)
    return true_ident_data, selected_data_id


def multi_parameter_box_plot(figure, outer_grid_object, number_parameters, parameter_value_info,
                             data_set_indices_to_plot=[]):
    """loop through multiple parameter data sets whose ident values are to be plotted as box plots"""
    inner_grid = gridspec.GridSpecFromSubplotSpec(number_parameters, 1,
                                                  subplot_spec=outer_grid_object, wspace=0.1, hspace=0.1)
    for i_parameter_number, i_parameter_info in enumerate(parameter_value_info):
        ax = plt.Subplot(figure, inner_grid[i_parameter_number])
        # get only identifying samples for each data set
        plot_data, selected_data_id  = get_true_boolean_data(i_parameter_info["raw sample ident data"],
                                          np.array(i_parameter_info["raw sample boolean data"]),
                                          data_set_indices_to_plot)
        multi_box_plot(figure, ax, plot_data, selected_data_id)
    return None


def plot_parameter_values(parameter_values, noise=0, data_sets_to_plot=[]):
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
    figure = plt.figure(figsize=(6, 4))
    outer_grid = gridspec.GridSpec(number_of_rows, number_of_columns, wspace=0.2, hspace=0.2)
    for i_flux in range(0, number_of_fluxes):
        if noise:
            # plot distribution of parameter values for all each data set that can identify parameters.
            # Distribution is plotted for all noise samples for each data sets used
            number_parameters = len(parameter_values["processed"][i_flux])
            try:
                multi_parameter_box_plot(figure, outer_grid[i_flux],
                                         number_parameters, parameter_values["processed"][i_flux],
                                         data_set_indices_to_plot=data_sets_to_plot)
            except TypeError:
                multi_parameter_box_plot(figure, outer_grid,
                                         number_parameters, parameter_values["processed"][i_flux],
                                         data_set_indices_to_plot=data_sets_to_plot)
        else:
            all_parameter_info = []
            all_parameter_names = []
            all_parameter_true_value = []
            number_parameters = len(parameter_values["processed"][i_flux])
            for i_parameter_info in parameter_values["processed"][i_flux]:
                all_parameter_info.append(i_parameter_info["sample mean"])
                all_parameter_names.append(i_parameter_info["parameter name"])
                all_parameter_true_value.append(i_parameter_info["true value"])
            try:
                axis_object = multi_box_plot(figure, outer_grid[i_flux], all_parameter_info, noise=0)
            except TypeError:
                axis_object = multi_box_plot(figure, outer_grid, all_parameter_info, noise=0)
            axis_object.scatter(range(1, number_parameters + 1), all_parameter_true_value,
                                  color='r', marker='o')
            axis_object.set_xticklabels(all_parameter_names)
            # Remove top axes and right axes ticks
            axis_object.get_xaxis().tick_bottom()
            axis_object.get_yaxis().tick_left()
    figure.show()
    plt.show()

    return None


def plot_experiment_data_dist(exp_xss, exp_fss, experiment_choice=()):
    exp_xss = np.array(exp_xss[0])
    number_experiments, number_mets = exp_xss.shape
    exp_fss = np.array(exp_fss[0])
    _, number_fluxes = exp_fss.shape
    # take only data from chosen experiment id
    if experiment_choice:
        chosen_exp_xss = exp_xss[experiment_choice, :]
        chosen_exp_fss = exp_fss[experiment_choice, :]
    else:
        chosen_exp_xss = exp_xss
        chosen_exp_fss = exp_fss
    # collate data for each metabolite
    all_metabolite_info = []
    for i_met in range(0, number_mets):
        all_metabolite_info.append(chosen_exp_xss[:, i_met])
    # collect data for each flux
    all_flux_info = []
    for i_flux in range(0, number_fluxes):
        all_flux_info.append(chosen_exp_fss[:, i_flux])

    # plot concentration and flux distributions
    f, axarr = plt.subplots(2, 1, figsize=(6, 4), dpi=100, facecolor='w', edgecolor='k')
    # plot concentration distributions
    bp = axarr[0].boxplot(all_metabolite_info)
    for whiskers in bp["whiskers"]:
        whiskers.set(color='k', linewidth=2)
    for flier in bp['fliers']:
        flier.set(marker='o', color='r', alpha=0.5)
    # set axis labels
    metabolite_names = variable_name('metabolite', range(0, number_mets))
    axarr[0].set_xticklabels(metabolite_names)
    # Remove top axes and right axes ticks
    axarr[0].get_xaxis().tick_bottom()
    axarr[0].get_yaxis().tick_left()

    # plot flux distribution
    bp = axarr[1].boxplot(all_flux_info)
    for whiskers in bp["whiskers"]:
        whiskers.set(color='k', linewidth=2)
    for flier in bp['fliers']:
        flier.set(marker='o', color='r', alpha=0.5)
    # set axis labels
    flux_names = variable_name('flux', range(0, number_fluxes))
    axarr[1].set_xticklabels(flux_names)
    # Remove top axes and right axes ticks
    axarr[1].get_xaxis().tick_bottom()
    axarr[1].get_yaxis().tick_left()
    plt.show()

    return None


def plot_numerical_parameter_estimates(all_parameter_info, noise=0):
    """get box plot of parameters estimated through numerical optimization"""
    if noise:
        pass
    else:
        number_parameters = len(all_parameter_info["values"])
        fig = plt.figure(figsize=(6, 4))
        outer_grid = gridspec.GridSpec(1, 1, wspace=0.2, hspace=0.2)
        try:
            axis_object = multi_box_plot(fig, outer_grid[0], all_parameter_info["values"], noise=0)
            axis_object.scatter(range(1, number_parameters + 1), all_parameter_info["true value"],
                                color='r', marker='o')
        except TypeError:
            axis_object = multi_box_plot(fig, outer_grid, all_parameter_info["values"], noise=0)
            axis_object.scatter(range(1, number_parameters + 1), all_parameter_info["true value"],
                                color='r', marker='o')

        axis_object.set_xticklabels(all_parameter_info["names"])
    plt.show()
    fig.show()
    return None


def plot_hist(fig_object, grid_object, distribution_data, mark_value=[], parameter_name=[],
              sharex_axis_object=[], sharey_axis_object=[]):
    """create subplot with grid object and plot histogram"""
    axis_object = plt.Subplot(fig_object, grid_object)
    hist_object = axis_object.hist(distribution_data)
    if parameter_name:
        axis_object.set_title(parameter_name, fontsize=18)
    if mark_value:
        mark_y_value = np.arange(0, np.max(hist_object[0]) + 2)
        mark_x_value = np.repeat(mark_value, len(mark_y_value))
        axis_object.plot(mark_x_value, mark_y_value, color='black', linestyle='dashed', linewidth=2)
    if sharex_axis_object and sharey_axis_object:
        fig_object.add_subplot(axis_object, sharex=sharex_axis_object, sharey=sharey_axis_object)
        plt.setp(axis_object.get_xticklabels(), visible=False)
        plt.setp(axis_object.get_yticklabels(), visible=False)
    elif sharex_axis_object:
        fig_object.add_subplot(axis_object, sharex=sharex_axis_object)
        plt.setp(axis_object.get_xticklabels(), visible=False)
        plt.ylabel('Frequency', fontsize=18)
    elif sharey_axis_object:
        fig_object.add_subplot(axis_object, sharey=sharey_axis_object)
        plt.setp(axis_object.get_yticklabels(), visible=False)
    else:
        plt.ylabel('Frequency', fontsize=18)
        fig_object.add_subplot(axis_object)
    return axis_object


def plot_box(fig_object, grid_object, distribution_data, sharex_axis_object=[], sharey_axis_object=[],
             vert_option=False):
    """create subplot with grid object and get box plot of distributions"""
    axis_object = plt.Subplot(fig_object, grid_object)
    bp_object = axis_object.boxplot(distribution_data, vert=vert_option)
    for whiskers in bp_object["whiskers"]:
        whiskers.set(color='k', linewidth=2)
    for flier in bp_object['fliers']:
        flier.set(marker='o', color='r', alpha=0.5)
    # Remove top axes and right axes ticks
    axis_object.get_xaxis().tick_bottom()
    axis_object.get_yaxis().tick_left()
    if sharex_axis_object and sharey_axis_object:
        fig_object.add_subplot(axis_object, sharex=sharex_axis_object, sharey=sharey_axis_object)
        plt.setp(axis_object.get_xticklabels(), visible=False)
        plt.setp(axis_object.get_yticklabels(), visible=False)
    elif sharex_axis_object:
        fig_object.add_subplot(axis_object, sharex=sharex_axis_object)
        plt.setp(axis_object.get_xticklabels(), visible=False)
    elif sharey_axis_object:
        fig_object.add_subplot(axis_object, sharey=sharey_axis_object)
        plt.setp(axis_object.get_yticklabels(), visible=False)
    else:
        fig_object.add_subplot(axis_object)
    return axis_object


def plot_violin(fig_object, grid_object, distribution_data, sharex_axis_object=[], sharey_axis_object=[]):
    """create subplot with grid object and get violin plot of distributions"""
    axis_object = plt.Subplot(fig_object, grid_object)
    violin_object = axis_object.violinplot(distribution_data, showmeans=True, showmedians=True)
    # for whiskers in bp_object["whiskers"]:
    #     whiskers.set(color='k', linewidth=2)
    # for flier in bp_object['fliers']:
    #     flier.set(marker='o', color='r', alpha=0.5)
    # Remove top axes and right axes ticks
    axis_object.get_xaxis().tick_bottom()
    axis_object.get_yaxis().tick_left()
    if sharex_axis_object and sharey_axis_object:
        fig_object.add_subplot(axis_object, sharex=sharex_axis_object, sharey=sharey_axis_object)
        plt.setp(axis_object.get_xticklabels(), visible=False)
        plt.setp(axis_object.get_yticklabels(), visible=False)
    elif sharex_axis_object:
        fig_object.add_subplot(axis_object, sharex=sharex_axis_object)
        plt.setp(axis_object.get_xticklabels(), visible=False)
    elif sharey_axis_object:
        fig_object.add_subplot(axis_object, sharey=sharey_axis_object)
        plt.setp(axis_object.get_yticklabels(), visible=False)
    else:
        fig_object.add_subplot(axis_object)
    return axis_object


def set_hist_box_axis_limits(hist_axis, box_axis):
    """set axis limits for hist and box plot of parameter values obtained from identifiability analysis"""
    box_x_lim = []
    hist_y_lim = []
    for i_box_axis, i_hist_axis in zip(box_axis, hist_axis):
        # get box axis x-limit
        box_x_lim.append(i_box_axis.get_xlim())
        # get hist y-limit
        hist_y_lim.append(i_hist_axis.get_ylim())

    # get maximimum of the histogram limits
    max_hist_lim = []
    min_hist_lim = []
    for i_hist in hist_y_lim:
        max_hist_lim.append(max(i_hist))
        min_hist_lim.append(min(i_hist))
    max_hist_lim = max(max_hist_lim)
    # min_hist_lim = min(min_hist_lim)

    # get maximums of box plot limits
    max_box_lim = []
    min_box_lim = []
    for i_box in box_x_lim:
        max_box_lim.append(max(i_box))
        min_box_lim.append(min(i_box))

    # set x-axis, y-axis limits
    for i_plot, (i_hist_axis, i_box_axis) in enumerate(zip(hist_axis, box_axis)):
        i_box_axis.set_xlim([0, max_box_lim[i_plot]])
        i_hist_axis.set_xlim([0, max_box_lim[i_plot]])
        i_hist_axis.set_ylim([0, max_hist_lim])
    plt.show()
    return None


def identify_column_row_numbers(number_of_plots):
    """get number of rows and columns based on number of plots"""
    if number_of_plots >= 3:
        number_of_rows = 3
    else:
        number_of_rows = 1
    if number_of_plots % number_of_rows != 0:
        number_of_columns = (number_of_plots + 1) / number_of_rows
    else:
        number_of_columns = number_of_plots / number_of_rows
    return number_of_rows, number_of_columns


def plot_parameter_value_hist(parameter_value, noise=0):
    """plot distribution of identified parameter values as a histogram"""
    number_of_fluxes = len(parameter_value["processed"])
    number_of_rows, number_of_columns = identify_column_row_numbers(number_of_fluxes)

    figure = plt.figure(figsize=(6, 4))
    outer_grid = gridspec.GridSpec(number_of_rows, number_of_columns, wspace=0.2, hspace=0.2)
    if noise:
        # this section is not complete - may not work May 1 2018
        for i_flux in range(0, number_of_fluxes):
            number_parameters = parameter_value["processed"][i_flux]
            inner_grid = gridspec.GridSpecFromSubplotSpec(2, number_parameters,
                                                          subplot_spec=outer_grid[i_flux], wspace=0.05, hspace=0.2)
            hist_axis_object = []
            for i_parameter_number, i_parameter_info in enumerate(parameter_value["processed"][i_flux]):
                hist_axis_object = plot_hist(figure, inner_grid[0, i_parameter_number], i_parameter_info["sample mean"],
                                             sharey_axis_object=hist_axis_object)
    else:
        for i_flux in range(0, number_of_fluxes):
            number_parameters = len(parameter_value["raw"][i_flux][0])
            inner_grid = gridspec.GridSpecFromSubplotSpec(2, number_parameters,
                                                          subplot_spec=outer_grid[i_flux], wspace=0.1, hspace=0.2)
            box_axis = []
            all_box_axis = []
            hist_axis = []
            all_hist_axis = []
            for i_parameter_number, i_parameter_info in enumerate(parameter_value["raw"][i_flux][0]):
                # box plot shares y-axis with other box plot(s)
                box_axis = plot_box(figure, inner_grid[1, i_parameter_number], i_parameter_info["found values"],
                                    sharey_axis_object=box_axis)
                all_box_axis.append(box_axis)
                # histogram shares y-axis with other histogram(s) and x-axis with box plot
                hist_axis = plot_hist(figure, inner_grid[0, i_parameter_number], i_parameter_info["found values"],
                                      parameter_name=i_parameter_info["parameter name"], mark_value=i_parameter_info["true values"][0],
                                      sharey_axis_object=hist_axis, sharex_axis_object=box_axis)
                all_hist_axis.append(hist_axis)

            set_hist_box_axis_limits(all_hist_axis, all_box_axis)
    return None


def plot_scatter(x_data, y_data, fig_object, grid_object):
    axis_object = plt.Subplot(fig_object, grid_object)
    axis_object.plot(x_data, y_data, linestyle='None', marker='o', markersize=2)
    fig_object.add_subplot(axis_object)
    return axis_object


def plot_line(x_data, y_data, fig_object, grid_object):
    axis_object = plt.Subplot(fig_object, grid_object)
    axis_object.plot(x_data, y_data, linestyle='--', color='black', linewidth=1.0)
    fig_object.add_subplot(axis_object)
    return axis_object


def plot_all_ss_estimates(original_ss, estimated_ss):
    """plot distribution of steady state concentrations/fluxes obtained from multiple parameter value estimates"""
    # create repetition of original values same as predicted value
    all_sample_original_ss = []
    for i_sample_ss, i_sample_exp_ss in zip(estimated_ss, original_ss):
        original_y_ss = [[i_perturbation_info] * len(i_sample_ss[i_perturbation])
                         for i_perturbation, i_perturbation_info in enumerate(i_sample_exp_ss)]
        all_sample_original_ss.append(original_y_ss)

    # plot details
    number_experiments = len(original_ss[0])
    number_plots = original_ss[0][0].shape[0]
    number_of_rows, number_of_columns = identify_column_row_numbers(number_plots)
    figure = plt.figure(figsize=(8, 6))
    outer_grid = gridspec.GridSpec(number_of_rows, number_of_columns, wspace=0.2, hspace=0.2)
    violin_axis_object = []
    box_axis_object = []
    for i_plot in range(0, number_plots):
        # generate x and y data
        x_data = [np.array(i_x_ss)[:, i_plot] for i_x_ss in all_sample_original_ss[0]]
        y_data = [np.array(j_x_ss)[:, i_plot] for j_x_ss in estimated_ss[0]]
        inner_grid = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer_grid[i_plot], wspace=0.1, hspace=0.2)
        # plot scatter on left side
        scatter_axis_object = plot_scatter(x_data, y_data, figure, inner_grid[0])
        x_data_diagonal = [i_experiment_data[i_plot] for i_experiment_data in original_ss[0]]
        scatter_axis_object.plot(x_data_diagonal, x_data_diagonal, linestyle='--', color='black', linewidth=1.0)
        # plot violin at second position
        violin_axis_object = plot_violin(figure, inner_grid[1], y_data, sharex_axis_object=violin_axis_object)
        # plot histogram on last position
        box_axis_object = plot_box(figure, inner_grid[2], y_data, vert_option=True,
                                   sharey_axis_object=violin_axis_object,
                                   sharex_axis_object=box_axis_object)
    return None
