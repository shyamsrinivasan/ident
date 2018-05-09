from generate_expdata import initialize_to_ss
from generate_expdata import perturb_parameters
from plot_ident_results import plot_all_ss_estimates
import copy
import matplotlib.pyplot as plt
import numpy as np


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
    all_ss_estimate = []
    all_dyn_estimate = []
    number_of_estimates = len(estimated_parameter)
    for j_parameter_id, (initial_ss, j_parameter_estimate) in enumerate(zip(estimate_initial_ss, estimated_parameter)):
        print("\nSimulation for parameter set {} of {}....\n".format(j_parameter_id + 1, number_of_estimates))
        estimate_perturbation_ss, _ = perturb_parameters(initial_ss, parameter_perturbation,
                                                                            cvode_options,
                                                                            j_parameter_estimate,
                                                                            number_of_samples,
                                                                            noise=noise,
                                                                            kinetics=kinetics,
                                                                            dynamic_plot=0,
                                                                            noise_std=noise_std)
        all_ss_estimate.append(estimate_perturbation_ss)
        # all_dyn_estimate.append(estimated_dyn)
    return all_ss_estimate, all_dyn_estimate


def form_dict_one_data_set(original_parameter, data_set_info):
    """form parameter dictionary from parameter estimates from one data set"""
    data_set_parameter_value = copy.deepcopy(original_parameter)
    for j_keys in data_set_info.keys():
        if j_keys in data_set_parameter_value:
            data_set_parameter_value[j_keys] = data_set_info[j_keys]
    return data_set_parameter_value


def form_dict_one_sample(original_parameter, sample_info, target_data=[]):
    """form parameter dictionary for estimated parameters from each
    sample of experimental data that contains multiple data sets"""
    all_data_set_values = []
    for data_iterator, (j_data_set_id, j_data_set_info) in \
            enumerate(zip(sample_info["data_id"], sample_info["parameter_value"])):
        if target_data:
            if data_iterator in set(target_data):
                # get new parameter dict for j_data_set
                j_data_set_parameter_value = form_dict_one_data_set(original_parameter, j_data_set_info)
            else:
                j_data_set_parameter_value = []
            all_data_set_values.append(j_data_set_parameter_value)
        else:
            # get new parameter dict for j_data_set
            j_data_set_parameter_value = form_dict_one_data_set(original_parameter, j_data_set_info)
            all_data_set_values.append(j_data_set_parameter_value)

    if target_data:
        select_data_set_values = [i_data_set_value
                                  for data_iterator, i_data_set_value in enumerate(all_data_set_values)
                                  if data_iterator in set(target_data)]
    else:
        select_data_set_values = all_data_set_values
    return select_data_set_values


def form_parameter_dict(original_parameter, extracted_parameters, target_samples=[], target_data=[]):
    """form parameter dictionaries from extracted parameters suitable for oden simulation"""
    number_of_samples = len(extracted_parameters)
    if not target_samples:
        # work with all samples
        all_sample_values = []
        for i_sample, i_sample_info in enumerate(extracted_parameters):
            print("Parameters for sample {} of {}:\n".format(i_sample + 1, number_of_samples))
            select_data_set_values = form_dict_one_sample(original_parameter, i_sample_info, target_data)
            all_sample_values.append(select_data_set_values)
    else:
        # work only with targeted samples
        all_sample_values = []
        pass
    return all_sample_values


def run_all_parameters(y0, cvode_options, original_parameter, extracted_parameter,
                       noise=0, kinetics=2, noise_std=0.05, target_data=[]):
    """validate parameter values based on obtained initial ss using parameters
        estimated using identifiability analysis"""
    # get all parameter sets in extracted parameter and form parameter dictionaries suitable for simulation
    all_sample_ode_parameters = form_parameter_dict(original_parameter, extracted_parameter, target_data=target_data)
    # get initial steady state information for original parameter set
    original_ss, original_dyn = initialize_to_ss(y0, cvode_options, original_parameter,
                                                 noise=noise, kinetics=kinetics, noise_std=noise_std)
    # simulate system with each estimated set of parameter values
    number_of_samples = len(all_sample_ode_parameters)
    all_sample_ss = []
    all_sample_dyn = []
    for j_sample, j_sample_parameter in enumerate(all_sample_ode_parameters):
        print("Parameters for sample {} of {}:\n".format(j_sample+1, number_of_samples))
        estimate_ss, estimate_dyn = run_initial_ss_simulation(y0, cvode_options, j_sample_parameter,
                                                              noise=noise, kinetics=kinetics, noise_std=noise_std)
        all_sample_ss.append(estimate_ss)
        all_sample_dyn.append(estimate_dyn)
    return original_ss, all_sample_ss, original_dyn, all_sample_dyn


def run_all_parameter_perturbation(y0, cvode_options, original_parameter, extracted_parameter,
                                   noise=0, kinetics=2, noise_std=0.05, target_data=[]):
    """run perturbation analysis for all estimated parameter data sets based on
    initial and perturbed steady states"""
    # get all parameter sets in extracted parameter and form parameter dictionaries suitable for simulation
    all_sample_ode_parameters = form_parameter_dict(original_parameter, extracted_parameter, target_data=target_data)

    # all parameter perturbations
    parameter_perturbation = [{"ac": 0}, {"ac": 1}, {"ac": 4}, {"ac": 9}, {"ac": -.1}, {"ac": -.5},
                              {"k1cat": .1}, {"k1cat": .5}, {"k1cat": 1}, {"k1cat": -.1}, {"k1cat": -.5},
                              {"V3max": .1}, {"V3max": .5}, {"V3max": 1}, {"V3max": -.1}, {"V3max": -.5},
                              {"V2max": .1}, {"V2max": .5}, {"V2max": 1}, {"V2max": -.1}, {"V2max": -.5}]

    # simulate system with each estimated set of parameter values
    number_of_samples = len(all_sample_ode_parameters)
    all_sample_ss = []
    all_sample_perturbation_ss = []
    all_sample_dyn = []
    all_sample_perturbation_dyn = []
    for j_sample, j_sample_parameter in enumerate(all_sample_ode_parameters):
        print("Parameters for sample {} of {}:\n".format(j_sample + 1, number_of_samples))
        # get initial system steady state for all estimated parameter values
        estimate_ss, estimate_dyn = run_initial_ss_simulation(y0, cvode_options, j_sample_parameter,
                                                              noise=noise, kinetics=kinetics, noise_std=noise_std)
        all_sample_ss.append(estimate_ss)
        all_sample_dyn.append(estimate_dyn)

        # run all perturbations for each estimated parameter value
        estimate_perturbation_ss, estimate_perturbation_dyn = run_perturbation_ss_simulation(estimate_ss, cvode_options,
                                                                                             j_sample_parameter,
                                                                                             parameter_perturbation,
                                                                                             number_of_samples=1,
                                                                                             noise=noise,
                                                                                             kinetics=kinetics,
                                                                                             noise_std=noise_std)
        all_sample_perturbation_ss.append(estimate_perturbation_ss)
        all_sample_perturbation_dyn.append(estimate_perturbation_dyn)

    # combine initial and perturbation ss for all samples
    all_sample_all_ss = []
    for j_sample_initial_ss, j_sample_perturbation_ss in zip(all_sample_ss, all_sample_perturbation_ss):
        # get all perturbation concentrations for each parameter estimate
        all_perturbation_y_ss = [[i_perturbation_info["y"] for i_perturbation_info in i_estimate_perturbation_ss]
                                 for i_estimate_perturbation_ss in j_sample_perturbation_ss]
        # get all perturbation fluxes for each parameter estimate
        all_perturbation_f_ss = [[i_perturbation_info["flux"] for i_perturbation_info in i_estimate_perturbation_ss]
                                 for i_estimate_perturbation_ss in j_sample_perturbation_ss]
        # get all perturbation ss id for each parameter estimate
        all_perturbation_ss_id = [[i_perturbation_info["ssid"] for i_perturbation_info in i_estimate_perturbation_ss]
                                  for i_estimate_perturbation_ss in j_sample_perturbation_ss]

        # get initial ss for each parameter estimate
        # all_initial_y_ss = [i_estimate_initial_ss["y"] for i_estimate_initial_ss in j_sample_initial_ss]
        # all_initial_f_ss = [i_estimate_initial_ss["flux"] for i_estimate_initial_ss in j_sample_initial_ss]
        # all_initial_ss_id = [i_estimate_initial_ss["ssid"] for i_estimate_initial_ss in j_sample_initial_ss]

        # combine initial ss with perturbation ss
        # number_parameter_estimates = len(all_initial_y_ss)
        # for i_estimate in range(0, number_parameter_estimates):
        #     all_perturbation_y_ss[i_estimate].insert(0, all_initial_y_ss[i_estimate])
        #     all_perturbation_f_ss[i_estimate].insert(0, all_initial_f_ss[i_estimate])
        #     all_perturbation_ss_id[i_estimate].insert(0, all_initial_ss_id[i_estimate])

        all_sample_all_ss.append({"y": all_perturbation_y_ss,
                                  "flux": all_perturbation_f_ss,
                                  "ssid": all_perturbation_ss_id})

    # combine initial and perturbation dynamics for all samples
    all_sample_all_dyn = []
    # for j_sample_initial_dyn, j_sample_perturbation_dyn in zip(all_sample_dyn, all_sample_perturbation_dyn):
    #     pass

    return all_sample_all_ss, all_sample_all_dyn


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


def validate_model(y0, cvode_options, original_parameter, extracted_parameter, experimental_data, ss=1, dyn=0,
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
    # compare new steady state with original experimental steady state
    if ss:
        # collect all ss values
        all_ss = collate_ss_values(all_sample_ss, experimental_data)
        plot_all_ss_estimates(all_ss["exp_y"], all_ss["y"])
        plot_all_ss_estimates(all_ss["exp_flux"], all_ss["flux"])
        # plot_ss_values(original_ss, all_ss, concentration=0, flux=1)

    # compare new dynamic values with original experimental dynamic values
    if dyn:
        # collate all dyn values
        pass

    return None


def validate_perturbation_ss():
    """validate estimated parameter values based on ss values obtained for different perturbations"""
    # get initial steady state information for all estimated parameter sets
    # get perturbation steady state information for all estimated parameter sets
    return None
