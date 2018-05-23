import copy
import numpy as np
import itertools as it
import pandas as pd
from collections import defaultdict
from generate_expdata import initialize_to_ss
from generate_expdata import perturb_parameters
from names_strings import variable_name


def run_initial_ss_simulation(y0, cvode_options, estimated_parameter, noise=0, kinetics=2, noise_std=0.05):
    """run ode simulations for original and all estimate parameter sets - may take a while"""
    # get initial steady state information for all estimated parameter sets
    all_ss_estimate = []
    all_dyn_estimate = []
    number_of_estimates = len(estimated_parameter)
    for j_parameter_id, j_parameter_estimate in enumerate(estimated_parameter):
        print("\nSimulation for parameter set {} of {}....\n".format(j_parameter_id + 1, number_of_estimates))
        estimated_ss, estimated_dyn = initialize_to_ss(y0, cvode_options, j_parameter_estimate,
                                                       noise=noise, kinetics=kinetics, noise_std=noise_std)
        all_ss_estimate.append(estimated_ss)
        all_dyn_estimate.append(estimated_dyn)
    return all_ss_estimate, all_dyn_estimate


def run_perturbation_ss_simulation(estimate_initial_ss, cvode_options, estimated_parameter, parameter_perturbation,
                                   number_of_samples=1, noise=0, kinetics=2, noise_std=0.05):
    """run ode simulations for all estimated parameter sets - may take a while 21 * number_estimated_parameters"""
    all_dyn_estimate = []
    all_ss_data_dict = defaultdict(list)
    empty_dict = {}
    number_of_estimates = len(estimated_parameter)
    for j_parameter_id, (initial_ss, j_parameter_estimate) in enumerate(zip(estimate_initial_ss, estimated_parameter)):
        print("\nSimulation for parameter set {} of {}....\n".format(j_parameter_id + 1, number_of_estimates))
        print("Parameter values: V3max: {}, K3fdp: {}, K3pep: {}\n".format(j_parameter_estimate["V3max"],
                                                                           j_parameter_estimate["K3fdp"],
                                                                           j_parameter_estimate["K3pep"]))
        estimate_name = 'estimate_{}'.format(j_parameter_id)
        estimate_perturbation_ss, _ = perturb_parameters(initial_ss, parameter_perturbation, cvode_options,
                                                         j_parameter_estimate, number_of_samples, noise=noise,
                                                         kinetics=kinetics, dynamic_plot=0, noise_std=noise_std)
        estimate_perturbation_ss["estimate_id"] = [estimate_name for _ in estimate_perturbation_ss["pep"]]
        for key, value in it.chain(empty_dict.items(), estimate_perturbation_ss.items()):
            for i_value in value:
                all_ss_data_dict[key].append(i_value)
        # all_dyn_estimate.append(estimated_dyn)
    return all_ss_data_dict, all_dyn_estimate


def form_dict_one_data_set(original_parameter, data_set_info):
    """form parameter dictionary from parameter estimates from one data set"""
    data_set_parameter_value = copy.deepcopy(original_parameter)
    for j_keys in data_set_info.keys():
        if j_keys in data_set_parameter_value:
            data_set_parameter_value[j_keys] = np.array([data_set_info[j_keys]])
    return data_set_parameter_value


def form_dict_one_sample(original_parameter, all_sample_info, select_combination_pos, target_data=[]):
    """form parameter dictionary for estimated parameters from each
    sample of experimental data that contains multiple data sets"""
    # get all parameter names and values for each data set
    all_data_set_values = []
    for i_sample_data_set_pair in select_combination_pos:
        parameter_values = [i_p_value[i_sample_data_set_pair] for i_p_value in all_sample_info["values"]]
        i_data_set_parameter_value = form_dict_one_data_set(original_parameter,
                                                            dict(zip(all_sample_info["names"], parameter_values)))
        all_data_set_values.append(i_data_set_parameter_value)

    return all_data_set_values


def form_parameter_dict(original_parameter, extracted_parameters, target_samples=[], target_data_set=[],
                        target_data=[]):
    """form parameter dictionaries from extracted parameters suitable for oden simulation"""
    # get all available sample names
    all_sample_names = [i_data_id[0] for i_data_id in extracted_parameters["data_sets"][0]]
    all_unique_sample_names = np.unique(np.array(all_sample_names)).tolist()
    select_sample_names = []
    if target_samples:
        # work only with targeted samples
        select_sample_names = set(target_samples).intersection(set(all_unique_sample_names))

    # get all available data sets
    all_data_set_ids = [i_data_id[1] for i_data_id in extracted_parameters["data_sets"][0]]
    all_unique_data_set_ids = np.unique(np.array(all_data_set_ids)).tolist()
    select_data_set_ids = []
    if target_data_set:
        # work only with selected data set ids
        select_data_set_ids = set(target_data_set).intersection(set(all_unique_data_set_ids))

    if select_sample_names and select_data_set_ids:
        # get all possible combinations
        number_data_sets = len(select_data_set_ids)
        all_repeated_sample_names = [j_sample_id for i_sample_id in select_sample_names
                                     for j_sample_id in [i_sample_id] * number_data_sets]
        number_samples = len(select_sample_names)
        all_repeated_data_set_ids = [j_data_set_id for j_data_set_id in [select_data_set_ids] * number_samples]
        possible_select_combinations = zip(all_repeated_sample_names, all_repeated_data_set_ids)

        # pick only available combinations of sample and data id
        select_combination = [i_select_combinations for i_select_combinations in possible_select_combinations
                              if i_select_combinations in set(extracted_parameters["data_sets"][0])]
    elif target_data:
        select_combination = [i_select_combinations for i_combination, i_select_combinations in
                              enumerate(extracted_parameters["data_sets"][0]) if i_combination in target_data]
    else:
        select_combination = extracted_parameters["data_sets"][0]

    select_combination_pos = [extracted_parameters["data_sets"][0].index(j_select_combo) for j_select_combo in
                              select_combination]

    all_sample_values = form_dict_one_sample(original_parameter, extracted_parameters, select_combination_pos)
    return all_sample_values, select_combination


def run_all_parameter_perturbation(y0, cvode_options, original_parameter, extracted_parameter,
                                   noise=0, kinetics=2, noise_std=0.05, target_data=[], target_data_set=[], target_sample=[]):
    """run perturbation analysis for all estimated parameter data sets based on
    initial and perturbed steady states"""
    # get all parameter sets in extracted parameter and form parameter dictionaries suitable for simulation
    all_sample_ode_parameters, select_combinations = form_parameter_dict(original_parameter, extracted_parameter,
                                                                         target_data=target_data,
                                                                         target_data_set=target_data_set,
                                                                         target_samples=target_sample)

    # all parameter perturbations
    parameter_perturbation = [{"wt": 0}, {"ac": 1}, {"ac": 4}, {"ac": 9}, {"ac": -.1}, {"ac": -.5},
                              {"k1cat": .1}, {"k1cat": .5}, {"k1cat": 1}, {"k1cat": -.1}, {"k1cat": -.5},
                              {"V3max": .1}, {"V3max": .5}, {"V3max": 1}, {"V3max": -.1}, {"V3max": -.5},
                              {"V2max": .1}, {"V2max": .5}, {"V2max": 1}, {"V2max": -.1}, {"V2max": -.5}]

    # simulate system with each estimated set of parameter values
    estimate_ss, estimate_dyn = run_initial_ss_simulation(y0, cvode_options, all_sample_ode_parameters,
                                                          noise=noise, kinetics=kinetics, noise_std=noise_std)

    # run all perturbations for each estimated parameter value
    estimate_perturbation_ss, estimate_perturbation_dyn = run_perturbation_ss_simulation(estimate_ss, cvode_options,
                                                                                         all_sample_ode_parameters,
                                                                                         parameter_perturbation,
                                                                                         number_of_samples=1,
                                                                                         noise=noise,
                                                                                         kinetics=kinetics,
                                                                                         noise_std=noise_std)
    # add sample and data set ids to perturbations data
    # unique_estimate_ids = np.unique(np.array(estimate_perturbation_ss["estimate_id"])).tolist()
    unique_experiment_ids = np.unique(np.array(estimate_perturbation_ss["experiment_id"]))

    all_sample_ids = [i_sample_name for i_estimate in select_combinations
                      for i_sample_name in [i_estimate[0]] * len(unique_experiment_ids)]
    all_data_set_ids = [i_data_set_id for i_estimate in select_combinations
                        for i_data_set_id in [i_estimate[1]] * len(unique_experiment_ids)]
    estimate_perturbation_ss.update({"sample_name": all_sample_ids, "data_set_id": all_data_set_ids})

    # combine initial and perturbation dynamics for all samples
    all_sample_all_dyn = []
    # for j_sample_initial_dyn, j_sample_perturbation_dyn in zip(all_sample_dyn, all_sample_perturbation_dyn):
    #     pass

    return estimate_perturbation_ss, all_sample_all_dyn


def collate_ss_values(ss_values, exp_ss_values):
    """collect and collate ss values from different data sets/perturbations or models/parameter sets"""
    # number_samples = len(ss_values)
    all_sample_y = []
    all_sample_f = []
    all_sample_exp_y = []
    all_sample_exp_f = []
    for i_sample, i_sample_info in enumerate(ss_values):
        # number_parameter_estimates = len(i_sample_info["y"])
        number_experiments = len(i_sample_info["y"][0])
        all_perturbation_y_ss = [[i_estimate_info[i_perturbation] for i_estimate_info in i_sample_info["y"]]
                                 for i_perturbation in range(0, number_experiments)]
        all_perturbation_f_ss = [[i_estimate_info[i_perturbation] for i_estimate_info in i_sample_info["flux"]]
                                 for i_perturbation in range(0, number_experiments)]
        all_sample_y.append(all_perturbation_y_ss)
        all_sample_f.append(all_perturbation_f_ss)
        all_sample_exp_y.append(exp_ss_values["y"][i_sample])
        all_sample_exp_f.append(exp_ss_values["flux"][i_sample])
    all_sample_info = {'y': all_sample_y,
                       'flux': all_sample_f,
                       'exp_y': all_sample_exp_y,
                       'exp_flux': all_sample_exp_f}
    return all_sample_info


def validate_model(y0, cvode_options, original_parameter, extracted_parameter, save_file_name=[], ss=1, dyn=0,
                   noise=0, kinetics=2, noise_std=0.05, target_data=[]):
    """vcalculate initial steady state for estimate parameter value"""
    # get initial and perturbation steady state information for original and all estimated parameter sets
    all_sample_ss, all_sample_dyn = run_all_parameter_perturbation(y0, cvode_options,
                                                                   original_parameter,
                                                                   extracted_parameter,
                                                                   noise=noise,
                                                                   kinetics=kinetics,
                                                                   noise_std=noise_std,
                                                                   target_data=target_data)
    # convert ss info to multi index data frame for storage and retrieval
    df_index_tuples = [(i_value_estimate, i_sample_id, i_data_set_id, i_exp_id)
                       for i_value_estimate, i_sample_id, i_data_set_id, i_exp_id in
                       zip(all_sample_ss["estimate_id"], all_sample_ss["sample_name"], all_sample_ss["data_set_id"],
                           all_sample_ss["experiment_id"])]
    multi_index_labels = ['estimate_id', 'sample_name', 'data_set_id', 'experiment_id']
    index = pd.MultiIndex.from_tuples(df_index_tuples, names=multi_index_labels)

    del all_sample_ss["estimate_id"]
    del all_sample_ss["sample_name"]
    del all_sample_ss["data_set_id"]
    del all_sample_ss["experiment_id"]

    # create data frame
    all_ss_df = pd.DataFrame(all_sample_ss, index=index, columns=all_sample_ss.keys())

    # save data frame to csv file
    if save_file_name:
        all_ss_df.to_csv(save_file_name, index_label=multi_index_labels)

    # compare new steady state with original experimental steady state
    if ss:
        # collect all ss values
        # all_ss = collate_ss_values(all_sample_ss, experimental_data)
        # plot_all_ss_estimates(all_ss["exp_y"], all_ss["y"])
        # plot_all_ss_estimates(all_ss["exp_flux"], all_ss["flux"])
        # plot_ss_values(original_ss, all_ss, concentration=0, flux=1)
        pass

    # compare new dynamic values with original experimental dynamic values
    if dyn:
        # collate all dyn values
        pass

    return all_ss_df


def get_variable_info(validate_df, exp_df, variable_type='metabolite'):
    """return all extracted info from validation and experimental df for
    desired type of variable (concentration/flux)"""

    idx = pd.IndexSlice

    # get all variable info from validation df
    # get all concentration names
    var_names = variable_name(variable_type)
    relevant_df = validate_df[var_names]
    # get all concentration values
    all_var_values = relevant_df.values
    all_var_list = [all_var_values[:, i_var] for i_var, _ in enumerate(var_names)]
    # scatter plot of concentration

    # number of parameter estiates
    number_estimates = len(relevant_df.index.levels[0].unique())

    # get all unique experiment info
    experiment_ids = relevant_df.index.levels[3].unique()
    # collect concentration info from all sample/data sets/experiments for distribution plots (violin/box)
    all_exp_var_values = [[relevant_df.loc[idx[:, :, :, i_experiment], idx[i_var]].values.tolist()
                          for i_experiment in experiment_ids] for i_var in var_names]
    experiment_ids = experiment_ids.tolist()

    # get steady state concentrations from original experimental data
    relevant_exp_df = exp_df[var_names]
    all_var_exp_values = relevant_exp_df.values
    all_var_exp_list = [all_var_exp_values[:, i_var] for i_var, _ in enumerate(var_names)]

    var_info = {"names": var_names,
                "values": [list(i_var) for i_var in all_var_list],
                "experiment_id": experiment_ids,
                "experiment_id_dist": all_exp_var_values,
                "experiment_values": [[x_val for i_val in j_var for x_val in [i_val]*number_estimates]
                                      for j_var in all_var_exp_list],
                "variable_type": variable_type}
    return var_info


def process_validation_info(validate_df, exp_df):
    """process all validation information from stored df for plotting"""
    idx = pd.IndexSlice
    # lexicographic ordering of df indices
    validate_df.sort_index(level="estimate_id", inplace=True)
    validate_df.sort_index(level="sample_name", inplace=True)
    validate_df.sort_index(level="data_set_id", inplace=True)
    validate_df.sort_index(level="experiment_id", inplace=True)

    # lexicographic ordering of exp df indices
    exp_df.sort_index(level='sample_name', inplace=True)
    exp_df.sort_index(level='experiment_id', inplace=True)

    # get all concentration (validation and original experiment for plotting
    all_c_info = get_variable_info(validate_df, exp_df, variable_type='metabolite')

    # get all flux (validation and original experiment)
    all_f_info = get_variable_info(validate_df, exp_df, variable_type='flux')

    return all_c_info, all_f_info
