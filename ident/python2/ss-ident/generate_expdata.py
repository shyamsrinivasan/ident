from generate_noisy_data import generate_no_noise_data
from generate_noisy_data import generate_noisy_data
from generate_noisy_data import run_no_noise_parameter_perturbation
from generate_noisy_data import run_noisy_parameter_perturbation
from plot_ident_results import plot_dynamic_sim_concentrations
from names_strings import variable_name
# import numpy as np
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


def create_ss_dict(ss_info, variable_type, noise=0):
    """create dictionary of all ss values variables-wise for use in creating data frames"""
    if noise:
        number_variables = len(ss_info[0][0])
    else:
        number_variables = len(ss_info[0])
    variable_name_info = [variable_name(variable_type, j_variable) for j_variable in range(0, number_variables)]

    variable_value_info = []
    if noise:
        for j_variable in range(0, number_variables):
            j_variable_info = []
            for i_sample_id, i_sample_info in enumerate(ss_info):
                for i_experiment_info in i_sample_info:
                    j_variable_info.append(i_experiment_info[j_variable])
            variable_value_info.append(j_variable_info)

        sample_name_info = []
        for i_sample_id, i_sample_info in enumerate(ss_info):
            for _ in i_sample_info:
                sample_name_info.append('sample_{}'.format(i_sample_id))
    else:
        for j_variable in range(0, number_variables):
            j_variable_info = []
            for i_experiment_id, i_experiment_info in enumerate(ss_info):
                j_variable_info.append(i_experiment_info[j_variable])
            variable_value_info.append(j_variable_info)

        sample_name_info = []
        i_sample_id = 0
        for _ in ss_info:
            sample_name_info.append('sample_{}'.format(i_sample_id))

    return variable_name_info, variable_value_info, sample_name_info


def create_other_value_dict(other_info, number_of_samples, noise=0):
    """create dictionary of other values based on number of samples
    to create consistent dict for data frame creation"""
    all_sample_final_ss = []
    for _ in range(0, number_of_samples):
        for i_experiment_value in other_info:
            all_sample_final_ss.append(i_experiment_value)
    return all_sample_final_ss


def perturb_parameters(initial_ss, parameter_perturbations, cvode_options, ode_parameter_values,
                       number_of_samples=1, noise=0, kinetics=2, dynamic_plot=0, noise_std=0.05):
    """perform parameter perturbations from given initial ss"""

    if noise:
        perturbation_options = {'ode_parameters': ode_parameter_values, 'cvode_options': cvode_options}
        final_ss, noisy_dynamic, experiment_info, _, dynamic_info = \
            run_noisy_parameter_perturbation(parameter_perturbations, initial_ss["y"][0], perturbation_options,
                                             kinetics=kinetics, number_of_samples=number_of_samples,
                                             noise_std=noise_std)
    else:
        perturbation_options = {'ode_parameters': ode_parameter_values, 'cvode_options': cvode_options}
        final_ss, no_noise_dynamic, experiment_info = \
            run_no_noise_parameter_perturbation(parameter_perturbations, initial_ss["y"], perturbation_options,
                                                kinetics=kinetics, plot_arg=dynamic_plot)

    # convert final_ss to dictionary suitable for data frame creation
    concentration_name, concentration_value, sample_name_info = create_ss_dict([i_ss["y"] for i_ss in final_ss],
                                                                               variable_type='metabolite', noise=noise)
    flux_name, flux_value, _ = create_ss_dict([i_ss["flux"] for i_ss in final_ss], variable_type='flux', noise=noise)

    # convert experiment_info to dict consistent with concentration_value and flux_value
    for i_field_name in experiment_info:
        new_field_value = create_other_value_dict(experiment_info[i_field_name], number_of_samples=number_of_samples, noise=noise)
        experiment_info[i_field_name] = new_field_value

    experiment_info.update(zip(concentration_name, concentration_value))
    experiment_info.update(zip(flux_name, flux_value))
    experiment_info.update({"sample_name": sample_name_info})

    # prepare list of column names for dataframe
    dict_fields = experiment_info.keys()
    # experiment_info_df = pd.DataFrame(experiment_info, columns=dict_fields)

    return experiment_info, dict_fields


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
        perturbation_info, info_dict_keys = perturb_parameters(initial_ss[0], parameter_perturbation,
                                                               cvode_options, ode_parameter_values,
                                                               number_of_samples, noise=noise, kinetics=kinetics,
                                                               dynamic_plot=perturbation_plot, noise_std=noise_std)
    except KeyError:
        perturbation_info, info_dict_keys = perturb_parameters(initial_ss, parameter_perturbation,
                                                               cvode_options, ode_parameter_values,
                                                               number_of_samples, noise=noise, kinetics=kinetics,
                                                               dynamic_plot=perturbation_plot, noise_std=noise_std)

    # convert info to multi index data frame for storage and retrieval
    df_index_tuples = [(i_value_sample, i_value_exp) for i_value_sample, i_value_exp in
                       zip(perturbation_info["sample_name"], perturbation_info["experiment_id"])]
    multi_index_labels = ['sample_name', 'experiment_id']
    index = pd.MultiIndex.from_tuples(df_index_tuples, names=multi_index_labels)
    del perturbation_info["sample_name"]
    del perturbation_info["experiment_id"]

    # create data frame
    all_ss_df = pd.DataFrame(perturbation_info, index=index, columns=perturbation_info.keys())

    return all_ss_df, multi_index_labels


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