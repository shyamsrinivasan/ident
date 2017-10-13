import numpy as np


def add_noise(concentration, flux=np.empty(shape=[0,6])):
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


def add_noise_dynamic():
    """function to add noise to input dynamic data"""
    pass
