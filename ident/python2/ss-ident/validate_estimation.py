from generate_expdata import initialize_to_ss


def calculate_initial_ss(y0, cvode_options, estimated_parameter, noise=0, kinetics=2, noise_std=0.05):
    """vcalculate initial steady state for estimate parameter value"""
    # get initial steady state information for all estimated parameter sets
    initial_ss, initial_dyn = initialize_to_ss(y0, cvode_options, estimated_parameter,
                                               noise=noise, kinetics=kinetics, noise_std=noise_std)

    # compare new steady state with original experimental steady state

    # get experimental system steady state data without noise using Convenience Kinetics for v3 (kinetics = 2)
    # exp_xss, exp_fss, exp_ssid, perturbation_details = \
    #     generate_expdata(y0, cvode_options, ode_parameter_values, noise=0, kinetics=2, dynamic_plot=0,
    #                      perturbation_plot=0)
    return None


def run_initial_ss_sim(y0, cvode_options, original_parameter, estimated_parameter, noise=0, kinetics=2, noise_std=0.05):
    """run ode simulations for original and all estimate parameter sets - may take a while"""
    # get initial steady state information for original parameter set
    initial_ss, initial_dyn = initialize_to_ss(y0, cvode_options, original_parameter,
                                               noise=noise, kinetics=kinetics, noise_std=noise_std)
    # get initial steady state information for all estimated parameter sets
    for j_parameter_estimate in estimated_parameter:
        estimated_ss, estimated_dyn = initialize_to_ss(y0, cvode_options, j_parameter_estimate, noise=noise, kinetics=kinetics, noise_std=noise_std)
    return None


def form_parameter_dict(original_parameters, extracted_parameters, target_samples=[]):
    """form parameter dictionaries from extracted parameters suitable for oden simulation"""
    number_of_samples = len(extracted_parameters)
    if not target_samples:
        # work with all samples
        all_sample_values = []
        for i_sample, i_sample_info in enumerate(extracted_parameters):
            print("Simulation for sample {} of {}:\n".format(i_sample, number_of_samples))
            all_data_set_values = []
            for data_iterator, (j_data_set_id, j_data_set_info) in \
                    enumerate(zip(i_sample_info["data_id"], i_sample_info["parameter_value"])):
                j_data_set_parameter_value = original_parameters[:]
                for j_keys in j_data_set_info.keys():
                    if j_keys in j_data_set_parameter_value.keys():
                        j_data_set_parameter_value[j_keys] = j_data_set_info[j_keys]
                all_data_set_values.append(j_data_set_parameter_value)
            all_sample_values.append(all_data_set_values)
    else:
        # work only with targeted samples
        pass

    return None


def validate_initial_ss(y0, cvode_options, original_parameter, extracted_parameter, noise=0, kinetics=2, noise_std=0.05):
    """validate parameter values based on obtained initial ss using parameters
                estimated using identifiability analysis"""
    # get all parameter sets in extracted parameter and form parameter dictionaries suitable for simulation
    return None


def validate_perturbation_ss():
    """validate estimated parameter values based on ss values obtained for different perturbations"""
    # get initial steady state information for all estimated parameter sets
    # get perturbation steady state information for all estimated parameter sets
    return None
