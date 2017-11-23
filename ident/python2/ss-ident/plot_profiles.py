import numpy as np
import matplotlib.pyplot as plt
plt.ion()


def plot_dynamic_course(time, concentration_data=np.array([]), flux_data=np.array([]), type=1):
    """use subplots from matplotlib sharing x per column and y per row"""
    if type == 1:  # concentrations
        if concentration_data.size:
            plt.plot(time, concentration_data, color='r')
            plt.xlabel('Time')
            plt.ylabel('Concentrations')
    elif type == 2:  # fluxes
        if flux_data.size:
            plt.plot(time, flux_data, color='g')
            plt.xlabel('Time')
            plt.ylabel('Fluxes')
    elif type == 3:  # both concentrations and fluxes in subplots
        f, axx = plt.subplots(2, sharex='col')
        axx[0].plot(time, concentration_data)
        axx[0].set_title('Concentrations')
        if flux_data.size:
            axx[1].plot(time, flux_data)
            axx[1].set_title('Fluxes')
    else:
        print('Incorrect plot type')

    #plt.show()
    #plt.pause(0.001)
    return None


def plot_multiple_dynamics(data):
    """plot multiple datasets passed in as a dictionary of values for time, y and flux"""
    for set_number in data:
        try:
            y_data = set_number["y"]
        except KeyError:
            y_data = set_number[0]
        try:
            f_data = set_number["flux"]
        except KeyError:
            f_data = set_number[1]
        try:
            time = set_number["time"]
        except KeyError:
            time = set_number[2]
        plot_dynamic_course(time, y_data, f_data, 3)
        plt.pause(0.0001)
    plt.ioff()
    plt.show()
    return None

