def ident_parameter_name(parameter_id, flux_name=(), flux_choice_id=0):
    parameter_list = ['V1max', 'K1ac (no enz)', 'k1cat', 'K1ac (enz)',
                      'V2max', 'K2pep',
                      'V3max', 'K3fdp', 'K3pep',
                      'V3max', 'K3fdp', 'K3pep',
                      'vemax', 'Kefdp',
                      'vemax', 'Kefdp']
    alter_list = {"flux1": [['V1max', 'K1ac (no enz)', 'k1cat', 'K1ac (enz)'],
                            ['V1max', 'K1ac (no enz)'],
                            ['k1cat', 'K1ac (enz)'],
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
