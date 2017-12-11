import numpy as np
import matplotlib.pyplot as plt
from kotte_model import ident_parameter_name


def plot_efficient_data():
    return None

def plot_identifiable_parameter(max_parameter):
    ident_data = max_parameter["data"]
    number_parameters = len(ident_data)
    y_pos = np.arange(number_parameters)
    # x_data = ident_data

    plt.rcdefaults()
    fig, ax = plt.subplots()
    ax.barh(y_pos, ident_data, align='center', color='green', ecolor='black')

    parameter_id = [ident_parameter_name(i) for i in range(0, number_parameters)]

    ax.set_yticks(y_pos)
    ax.set_yticklabels(parameter_id)
    ax.invert_yaxis()
    ax.set_xlabel('Number of Identifying Data Sets')
    ax.set_title('Parameter Identifiability')

    #max_ident_value = max_parameter["maximum"]
    #max_ident_data_id = max_parameter["id"]
    # number of datasets identifying maximum parameters
    #unique_parameter_number_identified = np.unique(np.array())

    #number_max_parameter_data = len(max_ident_data_id)
    return None

