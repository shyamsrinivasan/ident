from generate_noisy_data import generate_noisy_data
from generate_noisy_data import run_noisy_parameter_perturbation


def initialize_to_ss(y0, cvode_options, ode_parameter_values):
    """initialize system to ss from any y0 with given options and system parameters"""
    # get initial noisy system steady state
    initial_options = (cvode_options, ode_parameter_values)
    noisy_initial_ss, _, _, _ = generate_noisy_data(y0, initial_options, 1)
    return noisy_initial_ss


def perturb_parameters(initial_ss, parameter_perturbations, cvode_options, ode_parameter_values):
    """perform parameter perturbations from given initial ss"""
    perturbation_options = {'ode_parameters': ode_parameter_values, 'cvode_options': cvode_options}
    noisy_ss, noisy_dynamic, perturbation_details, _, dynamic_info = \
        run_noisy_parameter_perturbation(parameter_perturbations, initial_ss["y"], perturbation_options)
    return noisy_ss, perturbation_details


def generate_expdata(y0, cvode_options, ode_parameter_values):
    """generate noisy experimental data for kotte network to test identifiability"""
    # default parameter values
    # cvode_options = ('Newton', 'Adams', 1e-10, 1e-10, 200)
    # ode_parameter_values = np.array([.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, .1])

    # get initial noisy system steady state
    initial_ss = initialize_to_ss(y0, cvode_options, ode_parameter_values)

    # plot all dynamic courses
    # plot_multiple_dynamics(noisy_dynamic)
    # plt.close("all")
    # plot_multiple_dynamics(dynamic_info)
    # plt.close("all")

    # all parameter perturbations
    parameter_perturbation = [(14, 0), (14, 4), (14, 9),
                              (11, .1), (11, .5), (11, 1), (11, -.1), (11, -.5),
                              (12, .1), (12, .5), (12, 1), (12, -.1), (12, -.5),
                              (13, .1), (13, .5), (13, 1), (13, -.1), (13, -.5)]
    noisy_ss, perturbation_details = perturb_parameters(initial_ss, parameter_perturbation, cvode_options, ode_parameter_values)

    noisy_exp_xss = []
    noisy_exp_fss = []
    noisy_exp_ssid = []
    for ss_values in noisy_ss:
        noisy_exp_xss.append(ss_values["y"])
        noisy_exp_fss.append(ss_values["flux"])
        noisy_exp_ssid.append(ss_values["ssid"])

    return noisy_exp_xss, noisy_exp_fss, noisy_exp_ssid, perturbation_details


if __name__ == "__main__":
    import numpy as np

    y0 = np.array([5, 1, 1])
    # default parameter values
    cvode_options = ('Newton', 'Adams', 1e-10, 1e-10, 200)
    ode_parameter_values = np.array([.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, .1])

    # get initial noisy system steady state
    initial_ss = initialize_to_ss(y0, cvode_options, ode_parameter_values)

    # all parameter perturbations
    parameter_perturbation = [(14, 0), (14, 4), (14, 9),
                              (11, .1), (11, .5), (11, 1), (11, -.1), (11, -.5),
                              (12, .1), (12, .5), (12, 1), (12, -.1), (12, -.5),
                              (13, .1), (13, .5), (13, 1), (13, -.1), (13, -.5)]
    noisy_ss, perturbation_details = perturb_parameters(initial_ss, parameter_perturbation, cvode_options,
                                                        ode_parameter_values)
    noisy_exp_xss = []
    noisy_exp_fss = []
    noisy_exp_ssid = []
    for ss_values in noisy_ss:
        noisy_exp_xss.append(ss_values["y"])
        noisy_exp_fss.append(ss_values["flux"])
        noisy_exp_ssid.append(ss_values["ssid"])