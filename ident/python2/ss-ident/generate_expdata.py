from generate_noisy_data import generate_no_noise_data
from generate_noisy_data import generate_noisy_data
from generate_noisy_data import run_no_noise_parameter_perturbation
from generate_noisy_data import run_noisy_parameter_perturbation
from plot_ident_results import plot_dynamic_sim_concentrations


def initialize_to_ss(y0, cvode_options, ode_parameter_values, noise=1, kinetics=1):
    """initialize system to ss from any y0 with given options and system parameters"""
    if noise:
        # get initial noisy system steady state
        initial_options = (cvode_options, ode_parameter_values)
        # use non-noisy initial ss (wt) data
        noisy_initial_ss, noisy_initial_dyn, _, _ = generate_noisy_data(y0, initial_options, kinetics=kinetics)
        return noisy_initial_ss, noisy_initial_dyn
    else:
        # get initial noisy system steady state
        initial_options = (cvode_options, ode_parameter_values)
        initial_ss, initial_dyn, = generate_no_noise_data(y0, initial_options, kinetics=kinetics)
        return initial_ss, initial_dyn


def perturb_parameters(initial_ss, parameter_perturbations, cvode_options, ode_parameter_values,
                       number_of_samples=1, noise=1, kinetics=1, dynamic_plot=0):
    """perform parameter perturbations from given initial ss"""

    if noise:
        perturbation_options = {'ode_parameters': ode_parameter_values, 'cvode_options': cvode_options}
        noisy_ss, noisy_dynamic, perturbation_details, _, dynamic_info = \
            run_noisy_parameter_perturbation(parameter_perturbations, initial_ss["y"], perturbation_options,
                                             number_of_samples)
        return noisy_ss, perturbation_details
    else:
        perturbation_options = {'ode_parameters': ode_parameter_values, 'cvode_options': cvode_options}
        no_noise_ss, no_noise_dynamic, perturbation_details = \
            run_no_noise_parameter_perturbation(parameter_perturbations, initial_ss["y"], perturbation_options,
                                                kinetics=kinetics, plot_arg=dynamic_plot)
        return no_noise_ss, perturbation_details


def generate_expdata(y0, cvode_options, ode_parameter_values, number_of_samples=1, noise=1, kinetics=1,
                     dynamic_plot=0, perturbation_plot=0):
    """generate noisy experimental data for kotte network to test identifiability"""

    # get initial noisy system steady state
    initial_ss, initial_dyn = initialize_to_ss(y0, cvode_options, ode_parameter_values, noise, kinetics=kinetics)

    # plot all dynamic courses
    if dynamic_plot:
        plot_dynamic_sim_concentrations(initial_dyn, noise=noise)

    # all parameter perturbations
    parameter_perturbation = [{"ac": 0}, {"ac": 1}, {"ac": 4}, {"ac": 9}, {"ac": -.1}, {"ac": -.5},
                              {"k1cat": .1}, {"k1cat": .5}, {"k1cat": 1}, {"k1cat": -.1}, {"k1cat": -.5},
                              {"V3max": .1}, {"V3max": .5}, {"V3max": 1}, {"V3max": -.1}, {"V3max": -.5},
                              {"V2max": .1}, {"V2max": .5}, {"V2max": 1}, {"V2max": -.1}, {"V2max": -.5}]
    try:
        perturbed_ss, perturbation_details = perturb_parameters(initial_ss[0], parameter_perturbation,
                                                                cvode_options, ode_parameter_values,
                                                                number_of_samples, noise=noise, kinetics=kinetics,
                                                                dynamic_plot=perturbation_plot)
    except KeyError:
        perturbed_ss, perturbation_details = perturb_parameters(initial_ss, parameter_perturbation,
                                                                cvode_options, ode_parameter_values,
                                                                number_of_samples, noise=noise, kinetics=kinetics,
                                                                dynamic_plot=perturbation_plot)

    all_sample_exp_xss = []
    all_sample_exp_fss = []
    all_sample_exp_ssid = []
    for i_sample in range(0, number_of_samples):
        exp_xss = []
        exp_fss = []
        exp_ssid = []
        for ss_values in perturbed_ss:
            try:
                exp_xss.append(ss_values["y"][:, i_sample])
                exp_fss.append(ss_values["flux"][:, i_sample])
                exp_ssid.append(ss_values["ssid"][i_sample])
            except IndexError:
                exp_xss.append(ss_values["y"][:, ])
                exp_fss.append(ss_values["flux"][:, ])
                exp_ssid.append(ss_values["ssid"])
        all_sample_exp_xss.append(exp_xss)
        all_sample_exp_fss.append(exp_fss)
        all_sample_exp_ssid.append(exp_ssid)

    return all_sample_exp_xss, all_sample_exp_fss, all_sample_exp_ssid, perturbation_details


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