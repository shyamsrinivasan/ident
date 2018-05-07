from generate_expdata import initialize_to_ss
import copy
import numpy as np


def run_simulation(y0, cvode_options, estimated_parameter, noise=0, kinetics=2, noise_std=0.05):
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


def form_parameter_dict(original_parameter, extracted_parameters, target_samples=[]):
    """form parameter dictionaries from extracted parameters suitable for oden simulation"""
    number_of_samples = len(extracted_parameters)
    if not target_samples:
        # work with all samples
        all_sample_values = []
        for i_sample, i_sample_info in enumerate(extracted_parameters):
            print("Parameters for sample {} of {}:\n".format(i_sample+1, number_of_samples))
            all_data_set_values = []
            for data_iterator, (j_data_set_id, j_data_set_info) in \
                    enumerate(zip(i_sample_info["data_id"], i_sample_info["parameter_value"])):
                j_data_set_parameter_value = copy.deepcopy(original_parameter)
                for j_keys in j_data_set_info.keys():
                    if j_keys in j_data_set_parameter_value:
                        j_data_set_parameter_value[j_keys] = j_data_set_info[j_keys]
                all_data_set_values.append(j_data_set_parameter_value)
            all_sample_values.append(all_data_set_values)
    else:
        # work only with targeted samples
        all_sample_values = []
        pass
    return all_sample_values


def run_all_parameters(y0, cvode_options, original_parameter, extracted_parameter, noise=0, kinetics=2, noise_std=0.05):
    """validate parameter values based on obtained initial ss using parameters
        estimated using identifiability analysis"""
    # get all parameter sets in extracted parameter and form parameter dictionaries suitable for simulation
    all_sample_ode_parameters = form_parameter_dict(original_parameter, extracted_parameter)
    # get initial steady state information for original parameter set
    original_ss, original_dyn = initialize_to_ss(y0, cvode_options, original_parameter,
                                                 noise=noise, kinetics=kinetics, noise_std=noise_std)
    # simulate system with each estimated set of parameter values
    number_of_samples = len(all_sample_ode_parameters)
    all_sample_ss = []
    all_sample_dyn = []
    for j_sample, j_sample_parameter in enumerate(all_sample_ode_parameters):
        print("Parameters for sample {} of {}:\n".format(j_sample+1, number_of_samples))
        estimate_ss, estimate_dyn = run_simulation(y0, cvode_options, j_sample_parameter,
                                                   noise=noise, kinetics=kinetics, noise_std=noise_std)
        all_sample_ss.append(estimate_ss)
        all_sample_dyn.append(estimate_dyn)
    return original_ss, all_sample_ss, original_dyn, all_sample_dyn


def collate_ss_values(ss_values):
    """collect and collate ss values from different data sets/perturbations or models/parameter sets"""
    number_samples = len(ss_values)
    all_sample_info = []
    for i_sample_info in ss_values:
        all_xss = [j_ss_values["y"] for j_ss_values in i_sample_info]
        all_fss = [j_ss_values["flux"] for j_ss_values in i_sample_info]
        all_sample_info.append({'y': all_xss,
                                'flux': all_fss})
    return all_sample_info


def validate_model(y0, cvode_options, original_parameter, extracted_parameter, ss=1, dyn=0,
                   noise=0, kinetics=2, noise_std=0.05):
    """vcalculate initial steady state for estimate parameter value"""
    # get initial steady state information for original and all estimated parameter sets
    original_ss, all_sample_ss, original_dyn, all_sample_dyn = run_all_parameters(y0, cvode_options,
                                                                                  original_parameter,
                                                                                  extracted_parameter,
                                                                                  noise=noise,
                                                                                  kinetics=kinetics,
                                                                                  noise_std=noise_std)
    # compare new steady state with original experimental steady state
    if ss:
        # collect all ss values
        all_ss = collate_ss_values(all_sample_ss)
        for i_sample_ss in all_ss:
            original_y_ss = [original_ss["y"]] * len(i_sample_ss["y"])
            original_flux_ss = [original_ss["flux"]] * len(i_sample_ss["flux"])
            i_sample_ss["original_y_ss"] = original_y_ss
            i_sample_ss["original_flux_ss"] = original_flux_ss
        # comparison plots
        pass

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
