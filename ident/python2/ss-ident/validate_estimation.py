from generate_expdata import initialize_to_ss
import copy
import matplotlib.pyplot as plt


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
        estimate_ss, estimate_dyn = run_simulation(y0, cvode_options, j_sample_parameter,
                                                   noise=noise, kinetics=kinetics, noise_std=noise_std)
        all_sample_ss.append(estimate_ss)
        all_sample_dyn.append(estimate_dyn)
    return original_ss, all_sample_ss, original_dyn, all_sample_dyn


def run_all_parameter_perturbation():
    """run perturbation analysis for all estimated parameter data sets"""

    return None


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


def plot_scatter(x_data, y_data, axis_object):
    axis_object.scatter(x_data, y_data)
    return None


def plot_ss_values(original_value, estimated_value, concentration=1, flux=0):
    """plot steady state values of all concentrations and fluxes"""
    if concentration:
        # create repetition of original values same as predicted value
        all_sample_original_x_ss = []
        for i_sample_ss in estimated_value:
            original_y_ss = [original_value["y"]] * len(i_sample_ss["y"])
            all_sample_original_x_ss.append(original_y_ss)

        # plot details
        number_concentrations = original_value["y"].shape[0]
        f, ax = plt.subplots(number_concentrations, 1, figsize=(8, 6), dpi=100, facecolor='w',
                             edgecolor='k')
        for i_concentration, i_axis_object in enumerate(ax):
            x_data = [i_x_ss[i_concentration] for i_x_ss in all_sample_original_x_ss[0]]
            y_data = [j_x_ss[i_concentration] for j_x_ss in estimated_value[0]["y"]]
            plot_scatter(x_data, y_data, i_axis_object)

    if flux:
        all_sample_original_flux_ss = []
        for i_sample_ss in estimated_value:
            original_flux_ss = [original_value["flux"]] * len(i_sample_ss["flux"])
            all_sample_original_flux_ss.append(original_flux_ss)

        # plot details
        number_fluxes = original_value["flux"].shape[0]
        f_1, ax_1 = plt.subplots(number_fluxes, 1, figsize=(8, 6), dpi=100, facecolor='w',
                                 edgecolor='k')
        for i_flux, j_axis_object in enumerate(ax_1):
            x_data = [i_flux_ss[i_flux] for i_flux_ss in all_sample_original_flux_ss[0]]
            y_data = [j_flux_ss[i_flux] for j_flux_ss in estimated_value[0]["flux"]]
            plot_scatter(x_data, y_data, j_axis_object)

    return None


def validate_model(y0, cvode_options, original_parameter, extracted_parameter, ss=1, dyn=0,
                   noise=0, kinetics=2, noise_std=0.05, target_data=[]):
    """vcalculate initial steady state for estimate parameter value"""
    # get initial steady state information for original and all estimated parameter sets
    original_ss, all_sample_ss, original_dyn, all_sample_dyn = run_all_parameters(y0, cvode_options,
                                                                                  original_parameter,
                                                                                  extracted_parameter,
                                                                                  noise=noise,
                                                                                  kinetics=kinetics,
                                                                                  noise_std=noise_std,
                                                                                  target_data=target_data)
    # get perturbation steady state information for original and all estimated parameter sets


    # compare new steady state with original experimental steady state
    if ss:
        # collect all ss values
        all_ss = collate_ss_values(all_sample_ss)
        plot_ss_values(original_ss, all_ss, concentration=1)
        plot_ss_values(original_ss, all_ss, concentration=0, flux=1)

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
