import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from kotte_model import ident_parameter_name
from kotte_model import kotte_experiment_type_name
from kotte_model import kotte_variable_name
plt.ion()


def plot_on_axis_object(axis_obj, x_data, y_data, x_error, x_percent_mean, x_percent_std, x_max, noise=0):
    """given axis object plot given data on axis object along with all given annotations"""
    # plot bar graphs for x_data vs y_data with x_error error bars
    if noise:
        axis_obj.barh(y_data, x_data, xerr=x_error, align='center')
    else:
        # no error bars
        axis_obj.barh(y_data, x_data, align='center')
    # set x-axis limits
    axis_obj.set_xlim(0, x_max)
    # recompute data limits
    axis_obj.relim()
    axis_obj.autoscale_view()
    # set custom tick positions
    # axis_obj.xaxis.set_ticks(np.arange(0, x_max, 40))

    # annotate percentages onto each bar in the graph
    for j_bar, p in enumerate(axis_obj.patches):
        p.set_facecolor('blue')
        p.set_edgecolor('k')
        # set hatch pattern
        # p.set_hatch('/')

        # annotate percentages onto bar graphs
        if noise:
            y_annotation = "{:.2f} + {:.2f}%".format(x_percent_mean[j_bar],
                                                     x_percent_std[j_bar])
        else:
            y_annotation = "{:.2f}%".format(x_percent_mean[j_bar])
        axis_obj.annotate(y_annotation, (p.get_width(), p.get_y() + p.get_height()/2),
                          ha='center', va='center', xycoords='data')
    # set y axis ticks
    axis_obj.set_yticks(y_data)
    # set y axis tick labels (parameter names)
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
    plt.xticks(angles[:-1], x_data)

    # Draw ylabels
    axis_obj.set_rlabel_position(0)
    y_data_bool = [True if iy <= 0.0 else False for iy in y_data]
    # exception here when all y_data = 0.0
    if not all(y_data_bool):
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
    # exception here when all y_data = 0.0
    if not all(y_data_bool):
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
    metabolite_names = kotte_variable_name('metabolite', range(0, number_mets))
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
    flux_names = kotte_variable_name('flux', range(0, number_fluxes))
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
        figure = plt.figure(figsize=(6, 4))
        outer_grid = gridspec.GridSpec(1, 1, wspace=0.2, hspace=0.2)
        try:
            axis_object = multi_box_plot(figure, outer_grid[0], all_parameter_info["values"], noise=0)
            axis_object.scatter(range(1, number_parameters + 1), all_parameter_info["true value"],
                                color='r', marker='o')
        except TypeError:
            axis_object = multi_box_plot(figure, outer_grid, all_parameter_info["values"], noise=0)
            axis_object.scatter(range(1, number_parameters + 1), all_parameter_info["true value"],
                                color='r', marker='o')

        axis_object.set_xticklabels(all_parameter_info["names"])
    plt.show()
    return None


def plot_parameter_value_hist(parameter_value, noise=0):
    """plot distribution of identified parameter values as a histogram"""
    number_of_fluxes = len(parameter_value["processed"])
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
    if noise:
        # this section is not complete - may not work May 1 2018
        for i_flux in range(0, number_of_fluxes):
            number_parameters = parameter_value["processed"][i_flux]
            inner_grid = gridspec.GridSpecFromSubplotSpec(1, number_parameters,
                                                          subplot_spec=outer_grid[i_flux], wspace=0.1, hspace=0.1)
            for i_parameter_number, i_parameter_info in enumerate(parameter_value["processed"][i_flux]):
                ax = plt.Subplot(figure, inner_grid[i_parameter_number])
                ax.hist(i_parameter_info["sample mean"])
    else:
        for i_flux in range(0, number_of_fluxes):
            number_parameters = len(parameter_value["raw"][i_flux][0])
            inner_grid = gridspec.GridSpecFromSubplotSpec(1, number_parameters,
                                                          subplot_spec=outer_grid[i_flux], wspace=0.1, hspace=0.1)
            for i_parameter_number, i_parameter_info in enumerate(parameter_value["raw"][i_flux][0]):
                ax = plt.Subplot(figure, inner_grid[0, i_parameter_number])
                ax.hist(i_parameter_info["found values"])
    return None


def plot_all_initial_value_parameter_estimates():
    """plot distribution from multiple initial values"""
    return None



