import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.font_manager as fnt
import seaborn as sns
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


def plot_on_axis_object_box(axis_object, distribution_data, mark_value=[], vert_option=True):
    """get box plot of distributions"""
    # axis_object = plt.Subplot(fig_object, grid_object)
    bp_object = axis_object.boxplot(distribution_data, vert=vert_option)
    for whiskers in bp_object["whiskers"]:
        whiskers.set(color='k', linewidth=2)
    for flier in bp_object['fliers']:
        flier.set(marker='o', color='r', alpha=0.5)

    # plot actual value in red
    if mark_value:
        axis_object.scatter(range(1, len(distribution_data) + 1), mark_value,
                            color='r', marker='o')

    # Remove top axes and right axes ticks
    axis_object.get_xaxis().tick_bottom()
    axis_object.get_yaxis().tick_left()
    return None


def plot_on_axis_object_violin(axis_object, distribution_data):
    """create violin plot of distributions"""
    # axis_object = plt.Subplot(fig_object, grid_object)
    violin_object = axis_object.violinplot(distribution_data, showmeans=True, showmedians=True)

    # Remove top axes and right axes ticks
    axis_object.get_xaxis().tick_bottom()
    axis_object.get_yaxis().tick_left()
    return None


def plot_on_axis_object_hist(axis_object, distribution_data, mark_value=[], parameter_name=[], bins=[]):
    """create histogram of distribution data"""

    if bins:
        new_ax = sns.distplot(distribution_data, kde=False, ax=axis_object, bins=bins)
    else:
        new_ax = sns.distplot(distribution_data, kde=False, ax=axis_object)

    if parameter_name:
        axis_object.set_title(parameter_name, fontsize=14)

    if mark_value:
        axis_object.plot(np.repeat(mark_value, 2), new_ax.get_ylim(), color='black', linestyle='dashed', linewidth=0.8)
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


def set_hbar_axis_properties(axis_obj, y_data, y_tick_label, x_max, x_percent_mean, x_label=[], x_percent_std=[],
                             figure_title=[]):
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
            y_annotation = "{:.2f} + {:.2f}%".format(float(x_percent_mean[j_bar]),
                                                     float(x_percent_std[j_bar]))
        else:
            y_annotation = "{:.2f}%".format(float(x_percent_mean[j_bar]))
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


def set_polar_axis_limits(axis_object, y_limit):
    """set axis limits for polar plots"""
    axis_object.set_rmax(y_limit)
    y_ticks = range(0, int(y_limit) + 10, 10)
    axis_object.set_rticks(y_ticks)
    axis_object.set_rlabel_position(22.5)
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


def plot_scatter(axis_object, x_data, y_data):
    """plot scatter data"""
    axis_object.plot(x_data, y_data, linestyle='None', marker='o', markersize=2)
    return None


def set_scatter_axis_limits(axis_object, x_data, y_data):
    """set axis limits for scatter plots"""
    # set x-axis ticks
    if max(x_data) + 0.5 > 5:
        axis_object.set_xticks(np.arange(0, max(x_data) + 0.5, 1))
    else:
        axis_object.set_xticks(np.arange(0, max(x_data) + 0.5, 0.5))
    # set y-axis ticks
    if max(y_data) + 0.5 > 5:
        axis_object.set_yticks(np.arange(0, max(y_data) + 0.5, 1))
    else:
        axis_object.set_yticks(np.arange(0, max(y_data) + 0.5, 0.5))
    return None


def identifiability_plot(info_dict):
    """plot identifiability (number and percentage of data sets identifying each parameter)
    of every parameter in a given flux"""
    number_of_parameters = len(info_dict["names"])
    f, ax = plt.subplots(1, 1, figsize=(10, 8), dpi=100)
    x_data = info_dict["ident_mean"]
    y_data = range(0, number_of_parameters)
    gap = 0.05
    # y_add = [0, gap + 0.2, 2 * (gap + 0.2)]
    y_add = [i_y_data * (gap + 0.2) for i_y_data in y_data]
    x_error = info_dict["ident_std"]
    x_percent_mean = info_dict["ident_percent_mean"]
    x_percent_std = info_dict["ident_percent_std"]
    x_max = max(info_dict["total_data_sets"])
    x_label = ['Number of data combinations used for identification']
    plot_on_axis_object(ax, x_data, y_add, x_error=x_error)
    set_hbar_axis_properties(ax, y_add, y_tick_label=info_dict["names"], x_max=x_max, x_percent_mean=x_percent_mean,
                             x_percent_std=x_percent_std, x_label=x_label,
                             figure_title=info_dict["flux_name"][0]+' parameters')
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
                         figsize=(10, 8), dpi=100, gridspec_kw={"wspace": 0.2, "hspace": 0.2})
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
    for i_parameter, (_, i_parameter_name) in enumerate(zip(info_dict["exp_info"], info_dict["names"])):
        set_polar_axis_limits(ax[i_parameter], max(all_max_y_data))
        f_p = fnt.FontProperties(size=14, weight='demibold')
        ax[i_parameter].set_title(i_parameter_name, **{'font_properties': f_p})

    # add legend
    plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))
    return None


def parameter_values_plot(info_dict, original_values=(), violin=False, box=True, bins=[]):
    """plot distribution of parameter values as a box plot, violin plot and/or histogram"""
    number_parameters = len(info_dict["names"])
    if box:
        f1 = plt.figure(figsize=(10, 8), dpi=100, tight_layout=True)
        plot_grid = gridspec.GridSpec(2, number_parameters)

        # plot box plot
        box_axis = f1.add_subplot(plot_grid[0, :])
        if original_values:
            plot_on_axis_object_box(box_axis, info_dict["values"],
                                    mark_value=[original_values[i_name] for i_name in info_dict["names"]])
        else:
            plot_on_axis_object_box(box_axis, info_dict["values"],
                                    mark_value=[])
        box_axis.set_xticklabels(info_dict["names"])

        # plot histogram
        for i_parameter, (i_parameter_value, i_parameter_name) in enumerate(zip(info_dict["values"], info_dict["names"])):
            # parameter_name = info_dict["names"][i_parameter]
            hist_axis = f1.add_subplot(plot_grid[1, i_parameter])
            if original_values:
                plot_on_axis_object_hist(hist_axis, i_parameter_value, mark_value=original_values[i_parameter_name],
                                         parameter_name=i_parameter_name, bins=bins)
            else:
                plot_on_axis_object_hist(hist_axis, i_parameter_value, mark_value=[],
                                         parameter_name=i_parameter_name, bins=bins)

    if violin:
        f2 = plt.figure(figsize=(10, 8), dpi=100, tight_layout=True)
        plot_grid = gridspec.GridSpec(2, number_parameters)

        # plot box plot
        violin_axis = f2.add_subplot(plot_grid[0, :])
        plot_on_axis_object_violin(violin_axis, info_dict["values"])
        violin_axis.set_xticks(np.arange(1, len(info_dict["names"])+1))
        violin_axis.set_xticklabels(info_dict["names"])

        # plot histogram
        for i_parameter, (i_parameter_value, i_parameter_name) in enumerate(zip(info_dict["values"], info_dict["names"])):
            # parameter_name = info_dict["names"][i_parameter]
            hist_axis = f2.add_subplot(plot_grid[1, i_parameter])
            if original_values:
                plot_on_axis_object_hist(hist_axis, i_parameter_value, mark_value=original_values[i_parameter_name],
                                         parameter_name=i_parameter_name)
            else:
                plot_on_axis_object_hist(hist_axis, i_parameter_value, mark_value=[],
                                         parameter_name=i_parameter_name)
    return None


def validation_hist(values, names, figure_object, grid_objects):
    """plot histogram of validated variable distributions"""
    for i_variable, (i_var_value, i_var_name) in enumerate(zip(values, names)):
        hist_axis = figure_object.add_subplot(grid_objects[1, i_variable])
        plot_on_axis_object_hist(hist_axis, i_var_value, mark_value=[],
                                 parameter_name=i_var_name)
    return None


def validation_scatter(info_dict, grid_objects, figure_object):
    """plot validation values vs experimental values as scatter plots"""
    for i_variable, (_, _) in enumerate(zip(info_dict["values"], info_dict["names"])):
        scatter_axis = figure_object.add_subplot(grid_objects[2, i_variable])
        # scatter code
        plot_scatter(scatter_axis, info_dict["experiment_values"][i_variable],
                     info_dict["values"][i_variable])
        # line plot of experimental vs experimental
        scatter_axis.plot(info_dict["experiment_values"][i_variable],
                          info_dict["experiment_values"][i_variable],
                          **{'color': 'black', 'linestyle': 'dashdot', 'linewidth': 1.5})
        # set scatter axis ticks
        set_scatter_axis_limits(scatter_axis, x_data=info_dict["experiment_values"][i_variable],
                                y_data=info_dict["values"][i_variable])
    return None


def separate_validation_plot(info_dict, scatter=True, box=False, violin=True):
    """plot scatter, hist and box/violin plot for given variables in input dict"""
    number_variables = len(info_dict["names"])
    if scatter:
        f1 = plt.figure(figsize=(10, 8), dpi=100, tight_layout=True)
        plot_grid = gridspec.GridSpec(3, number_variables)
    else:
        f1 = plt.figure(figsize=(10, 8), dpi=100, tight_layout=True)
        plot_grid = gridspec.GridSpec(2, number_variables)

    # box plot
    if box:
        box_axis = f1.add_subplot(plot_grid[0, :])
        plot_on_axis_object_box(box_axis, info_dict["values"])
        box_axis.set_xticklabels(info_dict["names"])

        # plot histogram
        validation_hist(info_dict["values"], info_dict["names"], figure_object=f1, grid_objects=plot_grid)

        # scatter plot
        if scatter:
            validation_scatter(info_dict, figure_object=f1, grid_objects=plot_grid)

    # violin plot
    if violin:
        violin_axis = f1.add_subplot(plot_grid[0, :])
        plot_on_axis_object_violin(violin_axis, info_dict["values"])
        violin_axis.set_xticks(np.arange(1, len(info_dict["names"]) + 1))
        violin_axis.set_xticklabels(info_dict["names"])

        # plot histogram
        validation_hist(info_dict["values"], info_dict["names"], figure_object=f1, grid_objects=plot_grid)

        # scatter plot
        if scatter:
            validation_scatter(info_dict, figure_object=f1, grid_objects=plot_grid)
    return None


def experiment_based_validation(info_dict, box=False, violin=True, flux_id=()):
    """plot concentrations for different experiment separately to
    look at distribution within each experiments"""
    if flux_id:
        number_variables = len(flux_id)
    else:
        number_variables = len(info_dict["names"])
    f1 = plt.figure(figsize=(10, 8), dpi=100, tight_layout=True)
    plot_grid = gridspec.GridSpec(1, number_variables)
    if box:
        pass

    if violin:
        i_plot = 0
        for i_variable, i_var_name in enumerate(info_dict["names"]):
            if flux_id and (i_var_name in flux_id):
                violin_axis = f1.add_subplot(plot_grid[0, i_plot])
                plot_on_axis_object_violin(violin_axis, info_dict["experiment_id_dist"][i_variable])
                violin_axis.set_xticks(np.arange(1, len(info_dict["experiment_id"]) + 1))
                violin_axis.set_xticklabels(info_dict["experiment_id"])
                for tick in violin_axis.get_xticklabels():
                    tick.set_rotation(90)
                violin_axis.set_title(i_var_name)
                i_plot += 1
            elif not flux_id:
                violin_axis = f1.add_subplot(plot_grid[0, i_plot])
                plot_on_axis_object_violin(violin_axis, info_dict["experiment_id_dist"][i_variable])
                violin_axis.set_xticks(np.arange(1, len(info_dict["experiment_id"]) + 1))
                violin_axis.set_xticklabels(info_dict["experiment_id"])
                for tick in violin_axis.get_xticklabels():
                    tick.set_rotation(90)
                violin_axis.set_title(i_var_name)
                i_plot += 1
        plt.show()
    return None


def experiment_dist_plot(info_dict, box=False, violin=True):
    """plot distribution of experiment concentration and flux data"""
    number_variables = len(info_dict["names"])
    f1 = plt.figure(figsize=(10, 8), dpi=100, tight_layout=True)
    plot_grid = gridspec.GridSpec(2, number_variables)

    # box plot
    if box:
        box_axis = f1.add_subplot(plot_grid[0, :])
        plot_on_axis_object_box(box_axis, info_dict["experiment_values"])
        box_axis.set_xticklabels(info_dict["names"])

        # plot histogram
        validation_hist(info_dict["experiment_values"], info_dict["names"], figure_object=f1, grid_objects=plot_grid)

    # violin plot
    if violin:
        violin_axis = f1.add_subplot(plot_grid[0, :])
        plot_on_axis_object_violin(violin_axis, info_dict["experiment_values"])
        violin_axis.set_xticks(np.arange(1, len(info_dict["names"]) + 1))
        violin_axis.set_xticklabels(info_dict["names"])

        # plot histogram
        validation_hist(info_dict["experiment_values"], info_dict["names"], figure_object=f1, grid_objects=plot_grid)

    return None


def validation_plot(info_dict, concentration=True, flux=False, violin=True, box=False, flux_id=(), scatter=True,
                    experiment_dist=True):
    """plot values of concentrations and fluxes obtained from validation experiments"""
    if concentration:
        concentration_dict = info_dict["concentration"]
        # plot all concentrations together (irrespective of experiments)
        separate_validation_plot(concentration_dict, violin=violin, box=box, scatter=scatter)
        # plot experiment-wise
        experiment_based_validation(concentration_dict, violin=violin, box=box)
        # plot experiment only data distribution
        if experiment_dist:
            experiment_dist_plot(concentration_dict)

    if flux:
        flux_dict = info_dict["flux"]
        separate_validation_plot(flux_dict, violin=violin, box=box, scatter=scatter)
        # plot experiment-wise
        experiment_based_validation(flux_dict, violin=violin, box=box, flux_id=flux_id)
        # plot experiment only data distribution
        if experiment_dist:
            experiment_dist_plot(flux_dict)

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
