import numpy as np


def ident_parameter_name(parameter_id, flux_name=(), flux_choice_id=0):
    parameter_list = ['V1max', 'K1ac (no enz)', 'k1cat', 'K1ac (enz)',
                      'V2max', 'K2pep',
                      'V3max', 'K3fdp', 'K3pep',
                      'V3max', 'K3fdp', 'K3pep',
                      'vemax', 'Kefdp',
                      'vemax', 'Kefdp']
    alter_list = {"flux1": [['V1max', 'K1ac', 'k1cat', 'K1ac'],
                            ['V1max', 'K1ac'],
                            ['k1cat', 'K1ac'],
                            ['V1max', 'K1ac (no enz)', 'k1cat (exp 1)', 'k1cat (exp 2)']],
                  "flux2": [['V2max', 'K2pep']],
                  "flux3": [['V3max (1)', 'K3fdp (1)', 'K3pep (1)', 'V3max (2)', 'K3fdp (2)', 'K3pep (2)'],
                            ['V3max', 'K3fdp', 'K3pep'],
                            ['V3max', 'K3fdp', 'K3pep'],
                            ['V3max', 'K3fdp', 'K3pep'],
                            ['V3max', 'K3fdp', 'K3pep', 'L3fdp']],
                  "flux4": [],
                  "flux5": [['vemax (1)', 'Kefdp (1)', 'vemax (1)', 'Kefdp (2)'],
                            ['vemax', 'Kefdp'],
                            ['vemax', 'Kefdp']],
                  "flux6": [['V3max (1)', 'K3fdp', 'V3max (2)', 'K3pep'],
                            ['V3max', 'K3fdp'],
                            ['V3max', 'K3pep']]}
    if flux_name:
        try:
            parameter_name = [alter_list[name][choice_id][id]
                              for name, choice_id, id in zip(flux_name, flux_choice_id, parameter_id)]
        except TypeError:
            parameter_name = alter_list[flux_name][flux_choice_id][parameter_id]
    else:
        try:
            parameter_name = [parameter_list[id] for id in parameter_id]
        except TypeError:
            parameter_name = parameter_list[parameter_id]
    return parameter_name


def default_ident_parameter_name(flux_name, flux_choice):
    name_list = {"flux1": [['V1max', 'K1ac', 'k1cat', 'K1ac'],
                            ['V1max', 'K1ac'],
                            ['k1cat', 'K1ac'],
                            ['V1max', 'K1ac (no enz)', 'k1cat (exp 1)', 'k1cat (exp 2)']],
                  "flux2": [['V2max', 'K2pep']],
                  "flux3": [['V3max (1)', 'K3fdp (1)', 'K3pep (1)', 'V3max (2)', 'K3fdp (2)', 'K3pep (2)'],
                            ['V3max', 'K3fdp', 'K3pep'],
                            ['V3max', 'K3fdp', 'K3pep'],
                            ['V3max', 'K3fdp', 'K3pep'],
                            ['V3max', 'K3fdp', 'K3pep', 'L3fdp']],
                  "flux4": [],
                  "flux5": [['vemax (1)', 'Kefdp (1)', 'vemax (1)', 'Kefdp (2)'],
                            ['vemax', 'Kefdp'],
                            ['vemax', 'Kefdp']],
                  "flux6": [['V3max (1)', 'K3fdp', 'V3max (2)', 'K3pep'],
                            ['V3max', 'K3fdp'],
                            ['V3max', 'K3pep']]}
    try:
        select_parameter_name = [name_list[name][choice] for name, choice in zip(flux_name, flux_choice)]
    except TypeError:
        select_parameter_name = name_list[flux_name][flux_choice]
    return select_parameter_name


def parameter_name(parameter_id):
    parameter_list = ['K1ac', 'K3fdp', 'L3fdp', 'K3pep',
                      'K2pep', 'vemax', 'Kefdp', 'ne', 'd',
                      'V4max', 'k1cat', 'V3max', 'V2max', 'ac']
    try:
        return [parameter_list[id] for id in parameter_id]
    except TypeError:
        return parameter_list[parameter_id]


def true_parameter_values(flux_based=0, flux_name=(), flux_choice_id=0, parameter_id=()):
    if flux_based:
        alter_parameter_value = {"flux1": [{"K1ac": np.array([.1]), "k1cat": np.array([1]),
                                           "V1max": np.array([1])},
                                           {"V1max": np.array([1]), "K1ac": np.array([.1])},
                                           {"K1ac": np.array([.1]), "k1cat": np.array([1])},
                                           {"V1max": np.array([1]), "K1ac": np.array([.1]),
                                            "k1cat": np.array([1])}],
                                 "flux2": [{"K2pep": np.array([.3]), "V2max": np.array([1])}],
                                 "flux3": [{"K3fdp (1)": np.array([.1]), "K3pep (1)": np.array([.1]),
                                            "V3max (1)": np.array([1]), "K3fdp (2)": np.array([.1]),
                                            "K3pep (2)": np.array([.1]), "V3max (2)": np.array([1])},
                                           {"K3fdp": np.array([.1]), "K3pep": np.array([.1]),
                                            "V3max": np.array([1])},
                                           {"K3fdp": np.array([.1]), "K3pep": np.array([.1]),
                                            "V3max": np.array([1])},
                                           {"V3max": np.array([1]), "K3fdp": np.array([.1]), "K3pep": np.array([.1])},
                                           {"V3max": np.array([1]), "K3fdp": np.array([.1]), "K3pep": np.array([.1]),
                                            "L3fdp": np.array([4])}],
                                 "flux5": [{"vemax (1)": np.array([1.1]), "Kefdp (1)": np.array([.45]),
                                            "Kefdp (2)": np.array([.45])},
                                           {"vemax": np.array([1.1]), "Kefdp": np.array([.45])},
                                           {"vemax": np.array([1.1]), "Kefdp": np.array([.45])}],
                                 "flux6": [{"V3max (1)": np.array([1]), "K3fdp": np.array([.1]),
                                            "V3max (2)": np.array([1]), "K3pep": np.array([.1])},
                                           {"V3max": np.array([1]), "K3fdp": np.array([.1])},
                                           {"V3max": np.array([1]), "K3pep": np.array([.1])}]}
        try:
            parameter_value = [alter_parameter_value[name][choice_id][id]
                               for name, choice_id, id in zip(flux_name, flux_choice_id, parameter_id)]
        except TypeError:
            parameter_value = alter_parameter_value[flux_name][flux_choice_id][parameter_id]
        return parameter_value
    else:
        default_parameter_value = {"K1ac": np.array([.1]), "K3fdp": np.array([.1]), "L3fdp": np.array([4e6]),
                                   "K3pep": np.array([.1]), "K2pep": np.array([.3]), "vemax": np.array([1.1]),
                                   "Kefdp": np.array([.45]), "ne": np.array([2]), "d": np.array([.25]),
                                   "V4max": np.array([.2]), "k1cat": np.array([1]), "V3max": np.array([1]),
                                   "V2max": np.array([1]), "ac": np.array([.1])}
        return default_parameter_value


def kotte_experiment_type_name(experiment_id):
    experiment_type_name_list = ['wildtype acetate', 'acetate perturbation', 'k1cat perturbation',
                                 'V3max perturbation', 'V2max perturbation']
    try:
        return [experiment_type_name_list[index] for index in experiment_id]
    except TypeError:
        return experiment_type_name_list[experiment_id]


def experiment_name(experiment_id, experiment_details):
    try:
        parameter_changed = parameter_name([int(experiment_details["indices"][i, 0]) for i in experiment_id])
        parameter_value = [experiment_details["indices"][i, 1] for i in experiment_id]
    except TypeError:
        parameter_changed = parameter_name(int(experiment_details["indices"][experiment_id, 0]))
        parameter_value = experiment_details["indices"][experiment_id, 1]
    experiment_name_list = ['{} changes {}'.format(j_p_change, j_p_value)
                            for j_p_change, j_p_value in zip(parameter_changed, parameter_value)]
    return experiment_name_list


def variable_name(var_type, var_id=()):
    met_list = ['pep', 'fdp', 'E']
    flux_list = ['v1', 'v5', 'v3', 'v2', 'v4', 'v6']
    if var_type == 'metabolite':
        if hasattr(var_id, '__iter__'):
            if var_id:
                var_names = [met_list[j_var_id] for j_var_id in var_id]
            else:
                var_names = met_list
        else:
            var_names = met_list[var_id]

        # if var_id:
        #     try:
        #         var_names = [met_list[j_var_id] for j_var_id in var_id]
        #     except TypeError:
        #         var_names = met_list[var_id]
        # else:
        #     var_names = met_list
    elif var_type == 'flux':
        if hasattr(var_id, '__iter__'):
            if var_id:
                var_names = [flux_list[j_var_id] for j_var_id in var_id]
            else:
                var_names = flux_list
        else:
            var_names = flux_list[var_id]

        # if var_id:
        #     try:
        #         var_names = [flux_list[j_var_id] for j_var_id in var_id]
        #     except TypeError:
        #         var_names = flux_list[var_id]
        # else:
        #     var_names = flux_list
    else:
        var_names = []
    return var_names
