import numpy as np
from numpy.random import RandomState


def add_noise(data, default_shape, number_of_samples=1):
    """function to add noise to concentration or flux
    deafult_shape - tuple of deafult shape for concentration or flux in data
    number_of_samples - number of noisy samples requested"""

    # set random seed to get repeatable noisy data
    prng_1 = RandomState(12345678)
    rows, columns = data.shape
    if not rows or not columns:
        rows, columns = default_shape
    noise = prng_1.normal(0, .05, (rows, columns, number_of_samples))
    # concentration_noise = np.random.normal(0, .05, [rows, columns])
    noisy_data = np.zeros((rows, columns, number_of_samples))
    for i_sample in range(0, number_of_samples):
        noisy_data[:, :, i_sample] = data * (1 + noise[:, :, i_sample])
    return noisy_data


def add_noise_dynamic(concentration_dynamic=np.zeros(shape=[0,3]), flux_dynamic=np.zeros(shape=[0,6]),
                      number_of_samples=1):
    """function to add noise to input dynamic data.
    number_of_samples creates number_of_samples noisy data sets from original data set passed as input argument"""
    if concentration_dynamic.any():
        noisy_dynamic_concentration = add_noise(data=concentration_dynamic, default_shape=concentration_dynamic.shape,
                                                number_of_samples=number_of_samples)
    else:
        noisy_dynamic_concentration = concentration_dynamic

    if flux_dynamic.any():
        noisy_dynamic_flux = add_noise(data=flux_dynamic, default_shape=flux_dynamic.shape,
                                       number_of_samples=number_of_samples)
    else:
        noisy_dynamic_flux = flux_dynamic

    return noisy_dynamic_concentration, noisy_dynamic_flux
