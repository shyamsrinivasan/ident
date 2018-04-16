# generate noisy data for different types of kinetics used for regulated metabolic reactions
import matplotlib.pyplot as plt
from kotte_model import *
from simulate_ode import run_ode_sims
from add_noise import add_noise_dynamic
# from plot_profiles import plot_multiple_dynamics
from plot_ident_results import plot_dynamic_sims
from copy import deepcopy


def generate_data(y0, all_options, kinetics):
    _, ode_par_val = all_options
    if kinetics == 1:  # MWC kinetics
        time, y_dynamic = run_ode_sims(kotte_ode, y0, all_options, all_options[0][-1])[:2]
        # calculate dynamic flux data
        flux_dynamic = np.array(map(lambda x: kotte_flux(x, ode_par_val), y_dynamic))

    elif kinetics == 2:  # Convenience kinetics
        time, y_dynamic = run_ode_sims(kotte_ck_ode, y0, all_options, all_options[0][-1])[:2]
        # calculate dynamic flux
        flux_dynamic = np.array(map(lambda x: kotte_ck_flux(x, ode_par_val), y_dynamic))
    else:
        time = []
        y_dynamic = []
        flux_dynamic = []

    # info on bistability
    if y_dynamic[-1, 0] > y_dynamic[-1, 1]:
        bistable = 1
    elif y_dynamic[-1, 0] < y_dynamic[-1, 1]:
        bistable = 2
    else:
        bistable = 0

    # get ss info from dynamic data
    steady_state_info = {"y":y_dynamic[-1, :], "flux":flux_dynamic[-1, :], "ssid":bistable}
    dynamic_info = {"y":y_dynamic, "flux":flux_dynamic, "time":time, "ssid":bistable}

    return steady_state_info, dynamic_info


def run_parameter_perturbation(parameter_perturbation, y0, other_options):
    """run parameter perturbations based on tuple input parameter perturbation
    with first position of tuple being parameter id and second index being
    parameter value"""

    ode_parameters = other_options['ode_parameters']
    cvode_options = other_options['cvode_options']
    ss_data = []
    dynamic_data = []
    perturbed_parameter = []

    for index, p_value in enumerate(parameter_perturbation):
        print('Perturbation {}\n'.format(index + 1))
        parameter_id, parameter_change = p_value
        changed_ode_parameter = ode_parameters[:]
        changed_ode_parameter[parameter_id - 1] = changed_ode_parameter[parameter_id - 1] * (1 + parameter_change)
        all_options = (cvode_options, changed_ode_parameter)
        # generate data using MWC Kinetics
        ss_iter, dynamic_iter = generate_data(y0, all_options, 1)
        ss_data.append(ss_iter)
        dynamic_data.append(dynamic_iter)
        perturbed_parameter.append(changed_ode_parameter)

    return ss_data, dynamic_data, tuple(perturbed_parameter)


def generate_no_noise_data(y0, all_options, kinetics):
    steady_state_info, dynamic_info = generate_data(y0, all_options, kinetics=kinetics)
    # dynamic data
    concentration_dynamic = dynamic_info["y"]
    flux_dynamic = dynamic_info["flux"]
    # info on bistability
    if concentration_dynamic[-1, 0] > concentration_dynamic[-1, 1]:
        bistable = 1
    elif concentration_dynamic[-1, 0] < concentration_dynamic[-1, 1]:
        bistable = 2
    else:
        bistable = 0
    dynamic_info = {"y": concentration_dynamic, "flux": flux_dynamic,
                    "time": dynamic_info["time"], "ssid": bistable}
    steady_state_info = {"y": concentration_dynamic[-1, :],
                         "flux": flux_dynamic[-1, :], "ssid": bistable}
    return steady_state_info, dynamic_info


def generate_noisy_data(y0, all_options, kinetics, number_of_samples=1):
    steady_state_info, dynamic_info = generate_data(y0, all_options, kinetics)
    # dynamic data
    concentration_dynamic = dynamic_info["y"]
    flux_dynamic = dynamic_info["flux"]

    # add noise to dynamic data
    noisy_concentration_dynamic, noisy_flux_dynamic = \
        add_noise_dynamic(concentration_dynamic, flux_dynamic, number_of_samples)

    # get information for each sample separately
    # all_sample_noisy_dynamic_info = []
    # all_sample_noisy_steady_state_info = []
    all_sample_noisy_bistable = []
    for i_sample in range(0, number_of_samples):
        # info on bistability
        if noisy_concentration_dynamic[i_sample][-1, 0] > noisy_concentration_dynamic[i_sample][-1, 1]:
            bistable = 1
        elif noisy_concentration_dynamic[i_sample][-1, 0] < noisy_concentration_dynamic[i_sample][-1, 1]:
            bistable = 2
        else:
            bistable = 0
        all_sample_noisy_bistable.append(bistable)
    noisy_dynamic_info = {"y": noisy_concentration_dynamic, "flux": noisy_flux_dynamic,
                          "time": dynamic_info["time"], "ssid": all_sample_noisy_bistable}
    noisy_steady_state_info = {"y": [noisy_concentration_dynamic[i_sample][-1, :]
                                     for i_sample in range(0, number_of_samples)],
                               "flux": [noisy_flux_dynamic[i_sample][-1, :]
                                        for i_sample in range(0, number_of_samples)],
                               "ssid": all_sample_noisy_bistable}
        # all_sample_noisy_dynamic_info.append(noisy_dynamic_info)
        # all_sample_noisy_steady_state_info.append(noisy_steady_state_info)

    return noisy_steady_state_info, noisy_dynamic_info, \
           steady_state_info, dynamic_info


def run_no_noise_parameter_perturbation(parameter_perturbation, y0, other_options, plot_arg=0, kinetics=1):
    """run parameter perturbations based on tuple input parameter perturbation
        with first position of tuple being parameter id and second index being
        parameter value"""

    ode_parameters = other_options['ode_parameters']
    cvode_options = other_options['cvode_options']
    ss_info = []
    dynamic_info = []
    experiment_info_boolean = []
    # perturbed_parameter = np.zeros((len(parameter_perturbation), len(ode_parameters)))
    perturbed_parameter = []
    # perturbation_indices = np.zeros((len(parameter_perturbation), 2))
    perturbation_indices = []
    perturbation_ssid = np.zeros((len(parameter_perturbation), 2))

    for index, perturbation_value in enumerate(parameter_perturbation):
        print('Perturbation {}\n'.format(index + 1))
        parameter_name = perturbation_value.keys()[0]
        parameter_change = np.array(perturbation_value.values()[0])
        changed_ode_parameter = deepcopy(ode_parameters)
        changed_ode_parameter[parameter_name] = changed_ode_parameter[parameter_name] * (1 + parameter_change)
        all_options = (cvode_options, changed_ode_parameter)
        # generate data using MWC or Convinience Kinetics (choice specified using "kinetics" parameter
        ss_iter, dynamic_iter = generate_no_noise_data(y0, all_options, kinetics=kinetics)
        ss_info.append(ss_iter)
        dynamic_info.append(dynamic_iter)
        # initial ss info
        if y0[0] > y0[1]:
            perturbation_ssid[index, 0] = 1
        elif y0[0] < y0[1]:
            perturbation_ssid[index, 0] = 2
        # final ss info
        perturbation_ssid[index, 1] = ss_iter["ssid"]

        # perturbed_parameter[index, :] = changed_ode_parameter[:]
        perturbed_parameter.append(changed_ode_parameter)
        # perturbation_indices[index, :] = parameter_id - 1, parameter_change
        perturbation_indices.append(perturbation_value)
        # boolean_info = [False] * len(ode_parameters)
        # boolean_info[parameter_id - 1] = True
        # experiment_info_boolean.append(boolean_info)

    # plot all dynamic courses
    if plot_arg:
        plot_dynamic_sims(dynamic_info, multiple=1, concentrations=1, fluxes=1)
        # plt.close("all")

    perturbation_field_names = ['values', 'indices', 'original', 'boolean', 'ssid']
    perturbation_details = dict(zip(perturbation_field_names,
                                    [perturbed_parameter, perturbation_indices,
                                     ode_parameters, experiment_info_boolean,
                                     perturbation_ssid]))
    return ss_info, dynamic_info, perturbation_details


def run_noisy_parameter_perturbation(parameter_perturbation, y0, other_options, number_of_samples=1, plot_arg=0):
    """run parameter perturbations based on tuple input parameter perturbation
    with first position of tuple being parameter id and second index being
    parameter value"""

    ode_parameters = tuple(other_options['ode_parameters'])
    cvode_options = other_options['cvode_options']
    noisy_ss = []
    noisy_dynamic = []
    ss_info = []
    dynamic_info = []
    experiment_info_boolean = []
    perturbed_parameter = np.zeros((len(parameter_perturbation), len(ode_parameters)))
    perturbation_indices = np.zeros((len(parameter_perturbation), 2))
    perturbation_ssid = np.zeros((len(parameter_perturbation), 2))

    for index, p_value in enumerate(parameter_perturbation):
        print('Perturbation {}\n'.format(index+1))
        parameter_id, parameter_change = p_value
        changed_ode_parameter = np.array(ode_parameters[:])
        changed_ode_parameter[parameter_id-1] = changed_ode_parameter[parameter_id-1]*(1 + parameter_change)
        all_options = (cvode_options, changed_ode_parameter)
        # generate data using MWC Kinetics
        noisy_ss_iter, noisy_dynamic_iter, ss_iter, dynamic_iter = \
            generate_noisy_data(y0, all_options, 1, number_of_samples)

        noisy_ss.append(noisy_ss_iter)
        noisy_dynamic.append(noisy_dynamic_iter)
        ss_info.append(ss_iter)
        dynamic_info.append(dynamic_iter)
        # initial ss info
        if y0[0] > y0[1]:
            perturbation_ssid[index, 0] = 1
        elif y0[0] < y0[1]:
            perturbation_ssid[index, 0] = 2
        # final ss info
        perturbation_ssid[index, 1] = noisy_ss_iter["ssid"][0]

        perturbed_parameter[index, :] = changed_ode_parameter[:]
        perturbation_indices[index, :] = parameter_id-1, parameter_change
        boolean_info = [False]*len(ode_parameters)
        boolean_info[parameter_id-1] = True
        experiment_info_boolean.append(boolean_info)

    # plot all dynamic courses
    if plot_arg:
        plot_dynamic_sims(noisy_dynamic, multiple=1, concentrations=1, fluxes=1)
        plt.close("all")

    perturbation_field_names = ['values', 'indices', 'original', 'boolean', 'ssid']
    perturbation_details = dict(zip(perturbation_field_names,
                                    [perturbed_parameter, perturbation_indices,
                                     np.array(ode_parameters[:]), np.array(experiment_info_boolean[:]),
                                     perturbation_ssid]))
    return noisy_ss, noisy_dynamic, perturbation_details, ss_info, dynamic_info
