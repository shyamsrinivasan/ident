import numpy as np
import pandas as pd
import os.path
from simulate_data import arrange_experimental_data
from names_strings import true_parameter_values
from run_sims import ModelSim
import kotte_model


def create_experiment_data(save_file_name, noise=0, kinetics=2, number_samples=1, noise_std=0.05):
    """run generate_expdata and store resulting experimental data as csv file from data frame"""
    # generate no-noise experimental data for testing identifiability
    user_ode_opts = {'iter': 'Newton', 'discr': 'Adams', 'atol': 1e-10, 'rtol': 1e-10,
                     'time_points': 200, 'display_progress': True, 'verbosity': 30}
    # initial ss to begin all simulations from
    y0 = np.array([5, 1, 1])
    # get and set true parameter values, if available separately
    default_parameters = true_parameter_values()

    # create simulation object to simulate model with above parameters and initial conditions
    model_1 = ModelSim(kotte_model.kotte_ck_ode, kotte_model.kotte_ck_flux, noise=noise, **{'kinetics': kinetics,
                                                                                            'ode_opts': user_ode_opts,
                                                                                            't_final': 200,
                                                                                            'wt_y0': y0,
                                                                                            'i_parameter':
                                                                                                default_parameters,
                                                                                            'sample_size':
                                                                                                number_samples,
                                                                                            'noise_std': noise_std})
    # initial value determination for wt before perturbation
    wt_ss, wt_dynamic = model_1.run_initial_sim([default_parameters], ['default_parameters'])

    # all parameter perturbations
    parameter_perturbation = [{"wt": 0}, {"ac": 1}, {"ac": 4}, {"ac": 9}, {"ac": -.1}, {"ac": -.5},
                              {"k1cat": .1}, {"k1cat": .5}, {"k1cat": 1}, {"k1cat": -.1}, {"k1cat": -.5},
                              {"V3max": .1}, {"V3max": .5}, {"V3max": 1}, {"V3max": -.1}, {"V3max": -.5},
                              {"V2max": .1}, {"V2max": .5}, {"V2max": 1}, {"V2max": -.1}, {"V2max": -.5}]

    experiment_id = ['experiment_{}'.format(parameter_id) for parameter_id, _ in enumerate(parameter_perturbation)]
    experiment_details = model_1.change_parameter_values(parameter_perturbation)

    # call model.simulate to get initial (WT) steady state for all parameter sets strating from same y0
    model_1.sim_model(parameter=experiment_details, experiment_ids=experiment_id, initial_value=[wt_ss[0]['y']])

    # create dictionary suitable for writing to df
    experiment_df, multi_index_labels = model_1.create_df(parameter_perturbation, experiment_details)

    # get experimental system steady state data without noise using Convenience Kinetics for v3 (kinetics = 2)
    # experiment_df, multi_index_labels, dyn_df, dyn_labels = generate_expdata(y0, cvode_options, ode_parameter_values,
    #                                                                          noise=noise, kinetics=kinetics,
    #                                                                          dynamic_plot=0, perturbation_plot=0,
    #                                                                          number_of_samples=number_samples, noise_std=noise_std)

    # save data frame to csv file
    experiment_df.to_csv(save_file_name, index_label=multi_index_labels)
    # dyn_df.to_csv(save_dyn_file_name, index_label=dyn_labels)
    print(' Experiment Data written to file \n')

    return experiment_df


def retrieve_experimental_data_from_file(data_file_name, multi_index_label):
    """retrieve experimental data from csv file"""
    # read dataframe from csv file
    experiment_df = pd.read_csv(data_file_name, index_col=multi_index_label)
    return experiment_df


def create_data_for_analysis(experimental_data_df, experiment_choice_index, experiments_per_set, new_data_file_name):
    """arrange experimental data and create data frame and
    store for future use for use with identifiability analysis"""
    all_experiment_indices = ['experiment_0', 'experiment_1', 'experiment_2', 'experiment_3', 'experiment_4',
                              'experiment_5', 'experiment_6', 'experiment_7', 'experiment_8', 'experiment_9',
                              'experiment_10', 'experiment_11', 'experiment_12', 'experiment_13', 'experiment_14',
                              'experiment_15', 'experiment_16', 'experiment_17', 'experiment_18', 'experiment_19',
                              'experiment_20']

    experiment_index_choice = [all_experiment_indices[j_choice] for j_choice in experiment_choice_index]

    # get combinations of experimental datasets
    complete_data_dict, index_labels, experiment_choice = \
        arrange_experimental_data(experimental_data_df, experiments_per_set=experiments_per_set,
                                  experiment_choice=experiment_index_choice)
    # create data frame from data
    # same data frame to file
    df_tuples = [(j_sample, j_data_set, j_experiment) for j_sample, j_data_set, j_experiment in
                 zip(complete_data_dict["sample_name"], complete_data_dict["data_set_id"],
                     complete_data_dict["experiment_id"])]

    # create multi index
    index = pd.MultiIndex.from_tuples(df_tuples, names=index_labels)
    # remove redundant columns
    del complete_data_dict["sample_name"]
    del complete_data_dict["data_set_id"]
    del complete_data_dict["experiment_id"]

    # create data frame
    data_df = pd.DataFrame(complete_data_dict, index=index, columns=complete_data_dict.keys())

    # save data frame to csv file
    data_df.to_csv(new_data_file_name, index_label=index_labels)
    return data_df


def extract_and_create_data_for_analysis(original_data_file_name, original_index_labels, new_data_file_name,
                                         experiment_choice_index, experiments_per_set):
    """extract existing information from given file and create new file for future analysis"""
    experimental_df = retrieve_experimental_data_from_file(data_file_name=original_data_file_name,
                                                           multi_index_label=original_index_labels)
    # lexsort df
    experimental_df.sort_index(level=['sample_name', 'experiment_id'], inplace=True)
    create_data_for_analysis(experimental_df, experiment_choice_index, experiments_per_set, new_data_file_name)

    return None


def create_data_for_flux(flux_id, noise, number_samples=1, kinetics=2):
    """create and store data sets for each flux separately"""
    if flux_id == 'v1':
        experiment_choice_id = [0,
                                1, 2, 3, 4, 5,
                                11, 12, 13, 14, 15,
                                16, 17, 18, 19, 20]

        if noise:
            file_name = os.path.join(os.getcwd(), 'exp/experiments_noise_{}_samples'.format(number_samples))
            new_data_file_name = os.path.join(os.getcwd(),
                                              'exp/exp_v1_2_experiments_noise_{}_samples'.format(number_samples))
        else:
            if kinetics == 1:
                file_name = os.path.join(os.getcwd(), 'exp/experiments_mwc')
                new_data_file_name = os.path.join(os.getcwd(), 'exp/exp_v1_2_experiments_mwc')
            elif kinetics == 2:
                file_name = os.path.join(os.getcwd(), 'exp/experiments')
                new_data_file_name = os.path.join(os.getcwd(), 'exp/exp_v1_2_experiments')
        experiments_per_set = 2
    elif flux_id == 'v2':
        experiment_choice_id = [0,
                                1, 2, 3, 4, 5,
                                6, 7, 8, 9, 10,
                                11, 12, 13, 14, 15]

        if noise:
            file_name = os.path.join(os.getcwd(), 'exp/experiments_noise_{}_samples'.format(number_samples))
            new_data_file_name = os.path.join(os.getcwd(),
                                              'exp/exp_v2_2_experiments_noise_{}_samples'.format(number_samples))
        else:
            if kinetics == 1:
                file_name = os.path.join(os.getcwd(), 'exp/experiments_mwc')
                new_data_file_name = os.path.join(os.getcwd(), 'exp/exp_v2_2_experiments_mwc')
            elif kinetics == 2:
                file_name = os.path.join(os.getcwd(), 'exp/experiments')
                new_data_file_name = os.path.join(os.getcwd(), 'exp/exp_v2_2_experiments')
        experiments_per_set = 2
    elif flux_id == 'v3':
        experiment_choice_id = [0,
                                1, 2, 3, 4, 5,
                                6, 7, 8, 9, 10,
                                16, 17, 18, 19, 20]

        if noise:
            file_name = os.path.join(os.getcwd(), 'exp/experiments_noise_{}_samples'.format(number_samples))
            new_data_file_name = os.path.join(os.getcwd(),
                                              'exp/exp_v3_3_experiments_noise_{}_samples'.format(number_samples))
        else:
            if kinetics == 1:
                file_name = os.path.join(os.getcwd(), 'exp/experiments_mwc')
                new_data_file_name = os.path.join(os.getcwd(), 'exp/exp_v3_3_experiments_mwc')
            elif kinetics == 2:
                file_name = os.path.join(os.getcwd(), 'exp/experiments')
                new_data_file_name = os.path.join(os.getcwd(), 'exp/exp_v3_3_experiments')
        experiments_per_set = 3
    elif flux_id == 'v3a':
        experiment_choice_id = [0,
                                1, 2, 3, 4, 5,
                                6, 7, 8, 9, 10,
                                16, 17, 18, 19, 20]

        if noise:
            file_name = os.path.join(os.getcwd(), 'exp/experiments_noise_{}_samples'.format(number_samples))
            new_data_file_name = os.path.join(os.getcwd(),
                                              'exp/exp_v3_4_experiments_noise_{}_samples'.format(number_samples))
        else:
            if kinetics == 1:
                file_name = os.path.join(os.getcwd(), 'exp/experiments_mwc')
                new_data_file_name = os.path.join(os.getcwd(), 'exp/exp_v3_4_experiments_mwc')
            elif kinetics == 2:
                file_name = os.path.join(os.getcwd(), 'exp/experiments')
                new_data_file_name = os.path.join(os.getcwd(), 'exp/exp_v3_4_experiments')
        experiments_per_set = 4
    elif flux_id == 'v3b':
        experiment_choice_id = [0,
                                1, 2, 3, 4, 5,
                                6, 7, 8, 9, 10,
                                16, 17, 18, 19, 20]

        if noise:
            file_name = os.path.join(os.getcwd(), 'exp/experiments_noise_{}_samples'.format(number_samples))
            new_data_file_name = os.path.join(os.getcwd(),
                                              'exp/exp_v3_2_experiments_noise_{}_samples'.format(number_samples))
        else:
            if kinetics == 1:
                file_name = os.path.join(os.getcwd(), 'exp/experiments_mwc')
                new_data_file_name = os.path.join(os.getcwd(), 'exp/exp_v3_2_experiments_mwc')
            elif kinetics == 2:
                file_name = os.path.join(os.getcwd(), 'exp/experiments')
                new_data_file_name = os.path.join(os.getcwd(), 'exp/exp_v3_2_experiments')
        experiments_per_set = 2
    elif flux_id == 'v5':
        experiment_choice_id = [0,
                                1, 2, 3, 4, 5,
                                6, 7, 8, 9, 10,
                                11, 12, 13, 14, 15,
                                16, 17, 18, 19, 20]

        if noise:
            if number_samples == 5:
                file_name = os.path.join(os.getcwd(), 'exp/experiments_noise_5_samples')
                new_data_file_name = os.path.join(os.getcwd(), 'exp/exp_v5_2_experiments_noise_5_samples')
            elif number_samples == 500:
                file_name = os.path.join(os.getcwd(), 'exp/experiments_noise_500_samples')
                new_data_file_name = os.path.join(os.getcwd(), 'exp/exp_v5_2_experiments_noise_500_samples')
        else:
            if kinetics == 1:
                file_name = os.path.join(os.getcwd(), 'exp/experiments_mwc')
                new_data_file_name = os.path.join(os.getcwd(), 'exp/exp_v5_2_experiments_mwc')
            elif kinetics == 2:
                file_name = os.path.join(os.getcwd(), 'exp/experiments')
                new_data_file_name = os.path.join(os.getcwd(), 'exp/exp_v5_2_experiments')
        experiments_per_set = 2

    multi_index_labels = ['sample_name', 'experiment_id']
    extract_and_create_data_for_analysis(file_name, multi_index_labels, new_data_file_name,
                                         experiment_choice_id, experiments_per_set=experiments_per_set)
    return None


if __name__ == "__main__":
    file_name = os.path.join(os.getcwd(), 'exp/experiments')
    experiment_info_df = create_experiment_data(file_name, noise=0, kinetics=2)

    # file_name = os.path.join(os.getcwd(), 'exp/experiments_mwc')
    # experiment_info_df = create_experiment_data(file_name, noise=0, kinetics=1)

    # file_name = os.path.join(os.getcwd(), 'exp/experiments_noise_5_samples')
    # experiment_info_df = create_experiment_data(file_name, noise=1, kinetics=2, number_samples=5, noise_std=0.05)

    # file_name = os.path.join(os.getcwd(), 'exp/experiments_noise_500_samples')
    # experiment_info_df = create_experiment_data(file_name, noise=1, kinetics=2, number_samples=500, noise_std=0.05)
