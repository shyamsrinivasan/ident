import numpy as np
import matplotlib.pyplot as plt


def plot_dynamic_course(time, concentration_data=np.array([]), flux_data=np.array([]), type=1):
    """use subplots from matplotlib sharing x per column and y per row"""
    if type==1: # concentrations
        if concentration_data.size:
            plt.plot(time, ydatac)
            plt.xlabel('Time')
            plt.ylabel('Concentrations')
    elif type==2: # fluxes
        if flux_data.size:
            plt.plot(time, ydataf)
            plt.xlabel('Time')
            plt.ylabel('Fluxes')
    elif type==3: # both concentrations and fluxes in subplots
        f, axx = plt.subplots(2, sharex='col')
        axx[0].plot(time, ydatac)
        axx[0].set_title('Concentrations')
        if flux_data.size:
            axx[1].plot(time, ydataf)
            axx[1].set_title('Fluxes')
    else:
        print('Incorrect plot type')

    plt.show()
    return None

