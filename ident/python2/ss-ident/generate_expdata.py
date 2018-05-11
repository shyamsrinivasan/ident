from generate_noisy_data import generate_no_noise_data
from generate_noisy_data import generate_noisy_data
from generate_noisy_data import run_no_noise_parameter_perturbation
from generate_noisy_data import run_noisy_parameter_perturbation
from plot_ident_results import plot_dynamic_sim_concentrations
from kotte_model import kotte_variable_name
import pandas as pd


def initialize_to_ss(y0, cvode_options, ode_parameter_values, noise=0, kinetics=2, noise_std=0.05):
    """initialize system to ss from any y0 with given options and system parameters"""
    if noise:
        # get initial noisy system steady state
        initial_options = (cvode_options, ode_parameter_values)
        # use non-noisy initial ss (wt) data
        noisy_initial_ss, noisy_initial_dyn, _, _ = generate_noisy_data(y0, initial_options, kinetics=kinetics,
                                                                        noise_std=noise_std)
        return noisy_initial_ss, noisy_initial_dyn
    else:
        # get initial noisy system steady state
        initial_options = (cvode_options, ode_parameter_values)
        initial_ss, initial_dyn, = generate_no_noise_data(y0, initial_options, kinetics=kinetics)
        return initial_ss, initial_dyn


def perturb_parameters(initial_ss, parameter_perturbations, cvode_options, ode_parameter_values,
                       number_of_samples=1, noise=0, kinetics=2, dynamic_plot=0, noise_std=0.05):
    """perform parameter perturbations from given initial ss"""

    if noise:
        perturbation_options = {'ode_parameters': ode_parameter_values, 'cvode_options': cvode_options}
        noisy_ss, noisy_dynamic, perturbation_details, _, dynamic_info = \
            run_noisy_parameter_perturbation(parameter_perturbations, initial_ss["y"][0], perturbation_options,
                                             kinetics=kinetics, number_of_samples=number_of_samples,
                                             noise_std=noise_std)
        return noisy_ss, perturbation_details
    else:
        perturbation_options = {'ode_parameters': ode_parameter_values, 'cvode_options': cvode_options}
        no_noise_ss, no_noise_dynamic, perturbation_details = \
            run_no_noise_parameter_perturbation(parameter_perturbations, initial_ss["y"], perturbation_options,
                                                kinetics=kinetics, plot_arg=dynamic_plot)
        return no_noise_ss, perturbation_details


def create_variable_dict(ss_info, variable_type):
    """create dictionary of given variables for creating data frames for each experiment in each sample"""
    number_variables = len(ss_info[0][0])
    variable_name_info = [kotte_variable_name(variable_type, j_variable) for j_variable in range(0, number_variables)]

    all_variable_info = []
    for j_variable in range(0, number_variables):
        j_variable_info = []
        for i_sample_id, i_sample_info in enumerate(ss_info):
            for i_experiment_info in i_sample_info:
                j_variable_info.append(i_experiment_info[j_variable])
        all_variable_info.append(j_variable_info)

    sample_name_info = []
    for i_sample_id, i_sample_info in enumerate(ss_info):
        for _ in i_sample_info:
            sample_name_info.append('sample_{}'.format(i_sample_id))

    if variable_type == 'metabolite':
        variable_name_info.append('sample_name')
        all_variable_info.append(sample_name_info)
    return variable_name_info, all_variable_info


def generate_expdata(y0, cvode_options, ode_parameter_values, number_of_samples=1, noise=0, kinetics=2,
                     dynamic_plot=0, perturbation_plot=0, noise_std=0.05):
    """generate noisy experimental data for kotte network to test identifiability"""

    # get initial noisy system steady state
    initial_ss, initial_dyn = initialize_to_ss(y0, cvode_options, ode_parameter_values, noise,
                                               kinetics=kinetics, noise_std=noise_std)

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
                                                                dynamic_plot=perturbation_plot, noise_std=noise_std)
    except KeyError:
        perturbed_ss, perturbation_details = perturb_parameters(initial_ss, parameter_perturbation,
                                                                cvode_options, ode_parameter_values,
                                                                number_of_samples, noise=noise, kinetics=kinetics,
                                                                dynamic_plot=perturbation_plot, noise_std=noise_std)

    if noise:
        all_sample_exp_xss = [[i_perturbation_data["y"][i_sample] for i_perturbation_data in perturbed_ss]
                              for i_sample in range(0, number_of_samples)]
        all_sample_exp_fss = [[i_perturbation_data["flux"][i_sample] for i_perturbation_data in perturbed_ss]
                              for i_sample in range(0, number_of_samples)]
        all_sample_exp_ssid = [[i_perturbation_data["ssid"][i_sample] for i_perturbation_data in perturbed_ss]
                               for i_sample in range(0, number_of_samples)]
    else:
        all_sample_exp_xss = [[i_perturbation_data["y"] for i_perturbation_data in perturbed_ss]]
        all_sample_exp_fss = [[i_perturbation_data["flux"] for i_perturbation_data in perturbed_ss]]
        all_sample_exp_ssid = [[i_perturbation_data["ssid"] for i_perturbation_data in perturbed_ss]]

    all_sample_ss_info = {"y": all_sample_exp_xss, "flux": all_sample_exp_fss, "ssid": all_sample_exp_ssid}

    # convert data to pandas data frame for storage and retrieval later
    concentration_name, concentration_data = create_variable_dict(all_sample_ss_info["y"], variable_type='metabolite')
    flux_name, flux_data = create_variable_dict(all_sample_ss_info["flux"], variable_type='flux')

    all_variable_name = concentration_name[:]
    all_variable_info = concentration_data[:]
    for i_flux, i_flux_info in zip(flux_name, flux_data):
        all_variable_name.append(i_flux)
        all_variable_info.append(i_flux_info)
    all_ss_info = dict(zip(all_variable_name, all_variable_info))
    all_ss_df = pd.DataFrame(all_ss_info, columns=all_variable_name)
    all_ss_df.set_index(['sample_name', all_ss_df.index], inplace=True)

    return all_sample_ss_info, perturbation_details


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