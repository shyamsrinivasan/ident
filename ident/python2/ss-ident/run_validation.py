from run_sims import ModelSim
from names_strings import true_parameter_values
import numpy as np
import kotte_model


class ValidateSim(object):
    def __init__(self, rhs_fun, flux_fun, noise=0, **kwargs):
        self.rhs_fun = rhs_fun
        self.flux_fun = flux_fun
        self.noise = noise
        try:
            self.test_perturbations = kwargs['test_perturbations']
        except KeyError:
            self.test_perturbations = []
        # try:
        #     self.estimate_ids = kwargs['estimate_id']
        # except KeyError:
        #     self.estimate_id = 'estimate_x'
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
        self.wt_ss = []
        self.wt_dynamic = []
        self.perturbation_ss = []
        self.perturbation_dynamic = []
        self.estimated_parameters = []
        self.estimate_ids = []

    def create_parameter_list(self, estimate_info):
        """create name value pairs of estimated parameters followed by list of all parameters for use in validation"""
        # create dictionary (of length n_p) of parameters
        number_estimates = len(estimate_info['data_sets'])
        parameter_name_value_pair = [dict(zip(estimate_info['parameter_names'],
                                              [estimate_info['parameter_values'][i_parameter][i_estimate]
                                               for i_parameter, _ in enumerate(estimate_info['parameter_names'])]))
                                     for i_estimate in range(0, number_estimates)]

        # create list of all parameter values of size n_p with each of the above estimated values
        parameter_list = [self.i_parameter for _ in parameter_name_value_pair]
        for i_index, i_value in enumerate(parameter_list):
            for i_key in parameter_name_value_pair[i_index].keys():
                i_value[i_key] = parameter_name_value_pair[i_index][i_key]

        # create n_p simulation objects, each with one of the above n_p parameter vectors
        # valid_mod = []
        # for j_estimate in parameter_list:
        #     new_object = ValidateSim(kotte_model.kotte_ck_ode, kotte_model.kotte_ck_flux, **{'kinetics': kinetics,
        #                                                                                      'ode_opts': user_ode_opts,
        #                                                                                      't_final': 200,
        #                                                                                      'wt_y0': y0,
        #                                                                                      'i_parameter': j_estimate,
        #                                                                                      'sample_size': number_samples,
        #                                                                                      'noise_std': noise_std})
        #     valid_mod.append(new_object)
        estimate_ids = ['estimate_{}'.format(j_estimate) for j_estimate, _ in enumerate(parameter_list)]
        self.estimated_parameters = parameter_list
        self.estimate_ids = estimate_ids
        return self

    def validate_model(self):
        """run parallel validation method"""
        # run parallel initial sim (based on setup parallel ode for multiple parameter sets)
        # run parallel perturbation sim for each parameter estimate
        return None

if __name__ == '__main__':
    user_ode_opts = {'iter': 'Newton', 'discr': 'Adams', 'atol': 1e-10, 'rtol': 1e-10,
                     'time_points': 200, 'display_progress': True, 'verbosity': 30}
    # initial ss to begin all simulations from
    y0 = np.array([5, 1, 1])
    # get and set true parameter values, if available separately
    default_parameters = true_parameter_values()
    # all parameter perturbations
    parameter_perturbation = [{"wt": 0}, {"ac": 1}, {"ac": 4}, {"ac": 9}, {"ac": -.1}, {"ac": -.5},
                              {"k1cat": .1}, {"k1cat": .5}, {"k1cat": 1}, {"k1cat": -.1}, {"k1cat": -.5},
                              {"V3max": .1}, {"V3max": .5}, {"V3max": 1}, {"V3max": -.1}, {"V3max": -.5},
                              {"V2max": .1}, {"V2max": .5}, {"V2max": 1}, {"V2max": -.1}, {"V2max": -.5}]
    kinetics = 2
    number_samples = 1
    noise_std = 0.05
    validate_obj = ValidateSim(kotte_model.kotte_ck_ode, kotte_model.kotte_ck_flux, **{'kinetics': kinetics,
                                                                                       'ode_opts': user_ode_opts,
                                                                                       't_final': 200,
                                                                                       'wt_y0': y0,
                                                                                       'i_parameter':
                                                                                           default_parameters,
                                                                                       'test_perturbations':
                                                                                           parameter_perturbation,
                                                                                       'sample_size': number_samples,
                                                                                       'noise_std': noise_std})
    validate_obj.create_parameter_list(estimate_info)

