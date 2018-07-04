import kotte_model
from names_strings import true_parameter_values
from parallel_ode import setup_parallel_ode
from simulate_ode import setup_serial_ode
import numpy as np
from copy import deepcopy


class ModelSim(object):
    def __init__(self, rhs_fun, flux_fun, noise=0, **kwargs):
        self.rhs_fun = rhs_fun
        self.flux_fun = flux_fun
        self.noise = noise
        # kinetics
        try:
            self.kinetics = kwargs['kinetics']
        except KeyError:
            self.kinetics = 2
        # ode solver options
        try:
            self.ode_opts = kwargs['ode_opts']
        except KeyError:
            self.ode_opts = {'iter': 'Newton', 'discr': 'Adams', 'atol': 1e-10, 'rtol': 1e-10,
                             'time_points': 200, 'display_progress': True, 'verbosity': 30}
        # number of samples (applicable when noise=1)
        try:
            self.samples = kwargs['samples']
        except KeyError:
            self.sample = 1
        # simulation time horizon (for all sims)
        try:
            self.t_final = kwargs['t_final']
        except KeyError:
            self.t_final = 500
        # set wt initial value (to run all simulations)
        try:
            self.wt_y0 = kwargs['wt_y0']
        except KeyError:
            self.wt_y0 = []
        # default/initial parameter list
        try:
            self.i_parameter = kwargs['i_parameter']
        except KeyError:
            self.i_parameter = {}

        self.dynamic_info = []
        self.ss_info = []

    def run_initial_sim(self, parameter, parameter_ids=(), **kwargs):
        import pdb;pdb.set_trace()
        try:
            self.sim_model(parameter, parameter_ids, [kwargs['y0']])
        except KeyError:
            self.sim_model(self.i_parameter, parameter_ids, [self.wt_y0])
        return self

    def change_parameter_values(self, changed_parameter, default_parameter={}):
        """change default parameter by value metnioned in one of changed parameter"""
        if not default_parameter:
            default_parameter = self.i_parameter

        new_parameter_list = []
        for i_parameter_change in changed_parameter:
            new_parameter = deepcopy(default_parameter)
            parameter_name = list(i_parameter_change.keys())[0]
            parameter_change = np.array(list(i_parameter_change.values())[0])
            if parameter_name == 'wt':
                new_parameter['ac'] = new_parameter['ac'] * (1 + parameter_change)
            else:
                new_parameter[parameter_name] = new_parameter[parameter_name] * (1 + parameter_change)
            new_parameter_list.append(new_parameter)
        return new_parameter_list

    @staticmethod
    def collate_results(results, external_info, experiment_id, parameter_opt=1, i_value_opt=0):
        """collate results from all parallel simulations"""
        # import pdb; pdb.set_trace()

        collated_info = [{'info': i_value, 'id': j_value, 'y': results['y'][j_id],
                          'time': results['time'][j_id], 'flux': results['flux'][j_id]}
                         for i_value_id, i_value in zip(experiment_id, external_info)
                         for j_id, j_value in enumerate(results['id']) if j_value == i_value_id]

        if parameter_opt:
            experiment_info = {'parameter': [j_info['info'] for j_info in collated_info],
                               'time': [j_info['time'] for j_info in collated_info],
                               'y': [j_info['y'] for j_info in collated_info],
                               'id': [j_info['id'] for j_info in collated_info],
                               'flux': [j_info['flux'] for j_info in collated_info]}
        elif i_value_opt:
            experiment_info = {'initial_value': [j_info['info'] for j_info in collated_info],
                               'time': [j_info['time'] for j_info in collated_info],
                               'y': [j_info['y'] for j_info in collated_info],
                               'id': [j_info['id'] for j_info in collated_info],
                               'flux': [j_info['flux'] for j_info in collated_info]}
        else:
            experiment_info = {}
        return experiment_info

    def sim_model(self, parameter, experiment_ids, initial_value):
        """simulate model with defined rhs fun and given parameters, initial values and
        ode solver options (if given, else use default)"""

        # call parallel solver instance for multiple parameter/initial values
        import pdb;
        pdb.set_trace()
        if len(initial_value) > 1:
            sim_result = setup_parallel_ode(ode_rhs_fun=self.rhs_fun, flux_fun=self.flux_fun, parameters=parameter[0],
                                            y0=initial_value,
                                            t_final=self.t_final, experiment_id=experiment_ids, ode_opts=self.ode_opts,
                                            i_value_opt=1, parameter_opt=0)
            dynamic_info = sim_result
        elif len(parameter) > 1:
            sim_result = setup_parallel_ode(ode_rhs_fun=self.rhs_fun, flux_fun=self.flux_fun, parameters=parameter,
                                            y0=initial_value[0], t_final=self.t_final, experiment_id=experiment_ids,
                                            ode_opts=self.ode_opts,
                                            i_value_opt=0, parameter_opt=1)
            dynamic_info = self.collate_results(sim_result, parameter, experiment_ids)
        else:
            # use serial solver instance to solve for single parameter/initial values
            import pdb;pdb.set_trace()
            dynamic_info = setup_serial_ode(ode_fun=self.rhs_fun, y_initial=initial_value[0], t_final=self.t_final,
                                            opts=[self.ode_opts, parameter[0]])
            # collated_result = [{'info': parameter[0], 'id': 'experiment_0', 'y': dynamic_info['y'],
            #                     'time': dynamic_info['time']}]
            # calculate flux
            all_dyn_flux = []
            for i_y in dynamic_info['y']:
                all_dyn_flux.append(np.array(list(map(lambda x: self.flux_fun(x, parameter[0]), i_y))))
            dynamic_info['flux'] = all_dyn_flux

        # bistability info and ss values
        bistable = []
        ss_y = []
        ss_flux = []
        for j_experiment_y, j_experiment_flux in zip(dynamic_info['y'], dynamic_info['flux']):
            # info on bistability
            if j_experiment_y[-1, 0] > j_experiment_y[-1, 1]:
                bistable.append(1)
            elif j_experiment_y[-1, 0] < j_experiment_y[-1, 1]:
                bistable.append(2)
            else:
                bistable.append(0)
            # info on ss values
            ss_y.append(j_experiment_y[-1, :])
            ss_flux.append(j_experiment_flux[-1, :])
        self.dynamic_info = dynamic_info
        self.ss_info = {'y': ss_y, 'flux': ss_flux, 'ss_id': bistable}

        return self


if __name__ == '__main__':
    user_ode_opts = {'iter': 'Newton', 'discr': 'Adams', 'atol': 1e-10, 'rtol': 1e-10,
                     'time_points': 200, 'display_progress': True, 'verbosity': 30}
    # initial ss to begin all simulations from
    y0 = np.array([5, 1, 1])
    # get and set true parameter values, if available separately
    default_parameters = true_parameter_values()
    # create simulation object to simulate model with above parameters and initial conditions
    model_1 = ModelSim(kotte_model.kotte_ck_ode, kotte_model.kotte_ck_flux, noise=0, **{'kinetics': 2,
                                                                                        'ode_opts': user_ode_opts,
                                                                                        't_final': 200,
                                                                                        'wt_y0': y0,
                                                                                        'i_parameter':
                                                                                            default_parameters})
    # initial value determination for wt before perturbation
    model_1.run_initial_sim([default_parameters], 'default_parameters')

    # all parameter perturbations
    import pdb; pdb.set_trace()
    parameter_perturbation = [{"wt": 0}, {"ac": 1}, {"ac": 4}, {"ac": 9}, {"ac": -.1}, {"ac": -.5},
                              {"k1cat": .1}, {"k1cat": .5}, {"k1cat": 1}, {"k1cat": -.1}, {"k1cat": -.5},
                              {"V3max": .1}, {"V3max": .5}, {"V3max": 1}, {"V3max": -.1}, {"V3max": -.5},
                              {"V2max": .1}, {"V2max": .5}, {"V2max": 1}, {"V2max": -.1}, {"V2max": -.5}]

    # default set of parameters to begin simulations with
    # model_parameters = [{"K1ac": np.array([.1]), "K3fdp": np.array([.1]), "L3fdp": np.array([4e6]),
    #                      "K3pep": np.array([.1]), "K2pep": np.array([.3]), "vemax": np.array([1.1]),
    #                      "Kefdp": np.array([.45]), "ne": np.array([2]), "d": np.array([.25]), "V4max": np.array([.2]),
    #                      "k1cat": np.array([1]), "V3max": np.array([1]), "V2max": np.array([1]), "ac": np.array([.1])},
    #                     {"K1ac": np.array([.1]), "K3fdp": np.array([.1]), "L3fdp": np.array([4e6]),
    #                      "K3pep": np.array([.1]), "K2pep": np.array([.3]), "vemax": np.array([1.1]),
    #                      "Kefdp": np.array([.45]), "ne": np.array([2]), "d": np.array([.25]), "V4max": np.array([.2]),
    #                      "k1cat": np.array([1]), "V3max": np.array([1]), "V2max": np.array([1]), "ac": np.array([.01])},
    #                     {"K1ac": np.array([.1]), "K3fdp": np.array([.1]), "L3fdp": np.array([4e6]),
    #                      "K3pep": np.array([.1]), "K2pep": np.array([.3]), "vemax": np.array([1.1]),
    #                      "Kefdp": np.array([.45]), "ne": np.array([2]), "d": np.array([.25]), "V4max": np.array([.2]),
    #                      "k1cat": np.array([1]), "V3max": np.array([1]), "V2max": np.array([1]), "ac": np.array([.5])},
    #                     {"K1ac": np.array([.1]), "K3fdp": np.array([.1]), "L3fdp": np.array([4e6]),
    #                      "K3pep": np.array([.1]), "K2pep": np.array([.3]), "vemax": np.array([1.1]),
    #                      "Kefdp": np.array([.45]), "ne": np.array([2]), "d": np.array([.25]), "V4max": np.array([.2]),
    #                      "k1cat": np.array([1]), "V3max": np.array([1]), "V2max": np.array([1]), "ac": np.array([1])}]
    experiment_id = ['experiment_{}'.format(parameter_id) for parameter_id, _ in enumerate(parameter_perturbation)]
    experiment_details = model_1.change_parameter_values(parameter_perturbation)

    # call model.simulate to get initial (WT) steady state for all parameter sets strating from same y0
    initial_wt_result = model_1.run_initial_sim(parameter=experiment_details, parameter_ids=experiment_id)
    import pdb;pdb.set_trace()
    sim_result = model_1.sim_model(parameter=experiment_details, experiment_ids=experiment_id, initial_value=y0)

    # call model.simulate to get perturbed steady state for all
    # parameter perturbation from corresponding steady states y0
    # sim_perturbation_result = model_1.sim_model

    # call model.change_parameter
