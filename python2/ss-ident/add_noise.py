import numpy as np


def add_noise(concentration=np.empty(shape=[0,3]), flux=np.empty(shape=[0,6])):
    """function to add noise to input ss data"""
    if concentration:
        rows, columns = concentration.shape
        concentration_noise = np.random.normal(0, .01, [rows, columns])
        noisy_concentration = concentration * (1 + concentration_noise)
    else:
        noisy_concentration = concentration

    if flux:
        rows, columns = flux.shape
        flux_noise = np.random.normal(0, .01, [rows, columns])
        noisy_flux = flux * (1 + flux_noise)
    else:
        noisy_flux = flux

    return noisy_concentration, noisy_flux


def add_noise_dynamic(concentration_dynamic=np.zeros(shape=[0,3]), flux_dynamic=np.zeros(shape=[0,6])):
    """function to add noise to input dynamic data"""
    if concentration_dynamic.any():
        rows, columns = concentration_dynamic.shape
        concentration_noise = np.random.normal(0, .05, [rows, columns])
        noisy_dynamic_concentration = concentration_dynamic * (1 + concentration_noise)
    else:
        noisy_dynamic_concentration = concentration_dynamic

    if flux_dynamic.any():
        rows, columns = flux_dynamic.shape
        flux_noise = np.random.normal(0, .05, [rows, columns])
        noisy_dynamic_flux = flux_dynamic * (1 + flux_noise)
    else:
        noisy_dynamic_flux = flux_dynamic

    return noisy_dynamic_concentration, noisy_dynamic_flux
