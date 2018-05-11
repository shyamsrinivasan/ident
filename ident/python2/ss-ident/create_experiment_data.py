import numpy as np
import pandas as pd
from generate_expdata import generate_expdata


def create_experiment_data(file_name, noise=0, kinetics=2, number_samples=1, noise_std=0.05):
    """run generate_expdata and store resulting experimental data as csv file from data frame"""
    # generate no-noise experimental data for testing identifiability
    y0 = np.array([5, 1, 1])
    # default parameter values
    cvode_options = ('Newton', 'Adams', 1e-10, 1e-10, 200)
    ode_parameter_values = {"K1ac": np.array([.1]),
                            "K3fdp": np.array([.1]),
                            "L3fdp": np.array([4e6]),
                            "K3pep": np.array([.1]),
                            "K2pep": np.array([.3]),
                            "vemax": np.array([1.1]),
                            "Kefdp": np.array([.45]),
                            "ne": np.array([2]),
                            "d": np.array([.25]),
                            "V4max": np.array([.2]),
                            "k1cat": np.array([1]),
                            "V3max": np.array([1]),
                            "V2max": np.array([1]),
                            "ac": np.array([.1])}

    # get experimental system steady state data without noise using Convenience Kinetics for v3 (kinetics = 2)
    experiment_df, multi_index_labels = generate_expdata(y0, cvode_options, ode_parameter_values,
                                                         noise=noise, kinetics=kinetics,
                                                         dynamic_plot=0, perturbation_plot=0,
                                                         number_of_samples=number_samples, noise_std=noise_std)
    # save data frame to csv file
    experiment_df.to_csv(file_name, index_label=multi_index_labels)
    return experiment_df


def retrieve_experimental_data_from_file(file_name, multi_index_label):
    """retrieve experimental data from csv file"""
    # read dataframe from csv file
    experiment_df = pd.read_csv(file_name, index_col=multi_index_label)
    return experiment_df


def get_concentrations(data_frame, concentration_label):
    """retrieve all concentration from data frame"""
    row_indices = data_frame.index
    column_indices = data_frame.columns
    return None


def get_fluxes(data_frame, flux_label):
    """retrieve all fluxes from data frame"""
    return None


def retrieve_experimental_data(file_name, multi_index_lablel):
    """retrieve experimental data from csv file and collect concentrations, fluxes and
    other information from resulting data frame and pass as output arguments"""
    experimental_df = retrieve_experimental_data_from_file(file_name, multi_index_lablel)
    # get concentrations
    get_concentrations(experimental_df, ['pep', 'fdp', 'E'])
    # get fluxes
    # get experimental info
    return None


if __name__ == "__main__":
    file_name = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\ident\python2\ss-ident\experiments_noise_500_samples'
    experiment_info_df = create_experiment_data(file_name, noise=1, kinetics=2, number_samples=500, noise_std=0.05)
