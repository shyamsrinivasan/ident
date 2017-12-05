# generate noisy data for different types of kinetics used for regulated metabolic reactions
import matplotlib.pyplot as plt
from kotte_model import *
from simulate_ode import run_ode_sims
from add_noise import add_noise_dynamic
from plot_profiles import plot_multiple_dynamics
# from copy import deepcopy


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

    # get ss info from dynamic data
    steady_state_info = {"y":y_dynamic[-1, :], "flux":flux_dynamic[-1, :]}
    dynamic_info = {"y":y_dynamic, "flux":flux_dynamic, "time":time}

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


def generate_noisy_data(y0, all_options, kinetics):
    steady_state_info, dynamic_info = generate_data(y0, all_options, kinetics)
    # dynamic data
    concentration_dynamic = dynamic_info["y"]
    flux_dynamic = dynamic_info["flux"]

    # add noise to dynamic data
    noisy_concentration_dynamic, noisy_flux_dynamic = add_noise_dynamic(concentration_dynamic, flux_dynamic)
    noisy_dynamic_info = {"y":noisy_concentration_dynamic, "flux":noisy_flux_dynamic, "time":dynamic_info["time"]}
    noisy_steady_state_info = {"y":noisy_concentration_dynamic[-1, :], "flux":noisy_flux_dynamic[-1, :]}
    # noisy_concentration_steady_state =
    # noisy_flux_steady_state =

    return noisy_steady_state_info, noisy_dynamic_info, \
           steady_state_info, dynamic_info


def run_noisy_parameter_perturbation(parameter_perturbation, y0, other_options, plot_arg=0):
    """run parameter perturbations based on tuple input parameter perturbation
    with first position of tuple being parameter id and second index being
    parameter value"""

    ode_parameters = tuple(other_options['ode_parameters'])
    cvode_options = other_options['cvode_options']
    noisy_ss = []
    noisy_dynamic = []
    ss_info = []
    dynamic_info = []
    perturbed_parameter = np.zeros((len(parameter_perturbation), len(ode_parameters)))

    for index, p_value in enumerate(parameter_perturbation):
        print('Perturbation {}\n'.format(index+1))
        parameter_id, parameter_change = p_value
        changed_ode_parameter = np.array(ode_parameters[:])
        changed_ode_parameter[parameter_id-1] = changed_ode_parameter[parameter_id-1]*(1 + parameter_change)
        all_options = (cvode_options, changed_ode_parameter)
        # generate data using MWC Kinetics
        noisy_ss_iter, noisy_dynamic_iter, ss_iter, dynamic_iter = generate_noisy_data(y0, all_options, 1)
        noisy_ss.append(noisy_ss_iter)
        noisy_dynamic.append(noisy_dynamic_iter)
        ss_info.append(ss_iter)
        dynamic_info.append(dynamic_iter)
        perturbed_parameter[index, :] = changed_ode_parameter[:]

    # plot all dynamic courses
    if plot_arg:
        plot_multiple_dynamics(noisy_dynamic)
        plt.close("all")
        plot_multiple_dynamics(dynamic_info)
        plt.close("all")


    return noisy_ss, noisy_dynamic, perturbed_parameter, ss_info, dynamic_info
