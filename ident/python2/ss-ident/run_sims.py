import kotte_model
from names_strings import true_parameter_values
from parallel_ode import setup_parallel_ode
from simulate_ode import setup_serial_ode
import numpy as np
from copy import deepcopy
from add_noise import add_noise
from names_strings import variable_name


class ModelSim(object):
    def __init__(self, rhs_fun, flux_fun, noise=0, **kwargs):
        self.rhs_fun = rhs_fun
        self.flux_fun = flux_fun
        self.noise = noise
        try:
            self.sample_size = kwargs['sample_size']
        except KeyError:
            self.sample_size = 1

        try:
            self.noise_std = kwargs['noise_std']
        except KeyError:
            self.noise_std = 0.05

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

        self.wt_ss = []
        self.wt_dynamic = []
        self.dynamic_info = []
        self.ss_info = []
        self.noisy_dynamic_info = []
        self.noisy_ss_info = []

    def run_initial_sim(self, parameter, parameter_ids=(), **kwargs):
        try:
            self.sim_model(parameter, parameter_ids, [kwargs['y0']])
        except KeyError:
            self.sim_model(parameter, parameter_ids, [self.wt_y0])

        initial_ss = deepcopy(self.ss_info)
        self.wt_ss = initial_ss
        self.ss_info = []
        initial_dynamics = deepcopy(self.dynamic_info)
        self.wt_dynamic = initial_dynamics
        self.dynamic_info = []

        return initial_ss, initial_dynamics

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

        # if parameter_opt:
        #     experiment_info = {'parameter': [j_info['info'] for j_info in collated_info],
        #                        'time': [j_info['time'] for j_info in collated_info],
        #                        'y': [j_info['y'] for j_info in collated_info],
        #                        'id': [j_info['id'] for j_info in collated_info],
        #                        'flux': [j_info['flux'] for j_info in collated_info]}
        # elif i_value_opt:
        #     experiment_info = {'initial_value': [j_info['info'] for j_info in collated_info],
        #                        'time': [j_info['time'] for j_info in collated_info],
        #                        'y': [j_info['y'] for j_info in collated_info],
        #                        'id': [j_info['id'] for j_info in collated_info],
        #                        'flux': [j_info['flux'] for j_info in collated_info]}
        # else:
        #     experiment_info = {}
        return collated_info

    def sim_model(self, parameter, experiment_ids, initial_value):
        """simulate model with defined rhs fun and given parameters, initial values and
        ode solver options (if given, else use default)"""

        # call parallel solver instance for multiple parameter/initial values
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
            dynamic_info = setup_serial_ode(ode_fun=self.rhs_fun, y_initial=initial_value[0], t_final=self.t_final,
                                            opts=[self.ode_opts, parameter[0]])
            # collated_result = [{'info': parameter[0], 'id': 'experiment_0', 'y': dynamic_info['y'],
            #                     'time': dynamic_info['time']}]
            # calculate flux
            all_dyn_flux = []
            for index, i_experiment in enumerate(dynamic_info):
                dyn_flux = np.array(list(map(lambda x: self.flux_fun(x, parameter[0]), i_experiment['y'])))
                i_experiment.update({'flux': dyn_flux, 'id': experiment_ids[index]})

        # bistability info and ss values
        ss_data = []
        for j_experiment in dynamic_info:
            # info on bistability
            if j_experiment['y'][-1, 0] > j_experiment['y'][-1, 1]:
                stability_id = 1
            elif j_experiment['y'][-1, 0] < j_experiment['y'][-1, 1]:
                stability_id = 2
            else:
                stability_id = 0
            # info on ss values
            j_experiment.update({'ssid': stability_id})
            ss_data.append({'y': j_experiment['y'][-1, :], 'flux': j_experiment['flux'][-1, :], 'ssid': stability_id,
                            'id': j_experiment['id']})
        self.dynamic_info = dynamic_info
        self.ss_info = ss_data

        # add noise
        self.add_noise_dynamic()

        return self

    def add_noise_dynamic(self):
        """function to add noise to input dynamic data.
        number_of_samples creates number_of_samples noisy data sets from original data set passed as input argument"""

        if self.noise:
            all_noisy_dyn = []
            all_noisy_ss = []
            for i_value, i_ss_value in zip(self.dynamic_info, self.ss_info):
                noisy_y = add_noise(data=i_value['y'], default_shape=i_value['y'].shape,
                                    number_of_samples=self.sample_size, noise_std=self.noise_std)
                noisy_y_ss = add_noise(data=i_ss_value['y'], default_shape=i_ss_value['y'].shape,
                                       number_of_samples=self.sample_size, noise_std=self.noise_std)
                noisy_flux = add_noise(data=i_value['flux'], default_shape=i_value['flux'].shape,
                                       number_of_samples=self.sample_size, noise_std=self.noise_std)
                noisy_flux_ss = add_noise(data=i_ss_value['flux'], default_shape=i_ss_value['flux'].shape,
                                          number_of_samples=self.sample_size, noise_std=self.noise_std)
                all_noisy_dyn.append({'y': noisy_y, 'flux': noisy_flux})
                all_noisy_ss.append({'y': noisy_y_ss, 'flux': noisy_flux_ss})

            self.noisy_dynamic_info = all_noisy_dyn
            self.noisy_ss_info = all_noisy_ss
        return None

    def whatever(self, perturbations, experiment_details):
        """create dictionary from results suitable for writing to df"""
        # convert perturbation details to dictionary suitable for data_frame creation
        parameter_name = [list(i_perturbation_info.keys())[0] for i_perturbation_info in perturbations]
        parameter_change = [np.array(i_perturbation_info.values()[0])
                            for i_perturbation_info in perturbations]
        parameter_value = [np.array(i_parameter_value_dict[i_parameter_name][0]) if i_parameter_name != 'wt'
                           else np.array(i_parameter_value_dict['ac'][0])
                           for i_parameter_name, i_parameter_value_dict in zip(parameter_name, experiment_details)]
        essential_parameter_value = [np.array(i_perturbation_info["ac"][0]) for i_perturbation_info in
                                     experiment_details]
        parameter_change_percentage = [i_parameter_change * 100 for i_parameter_change in parameter_change]
        initial_value_ss_id = [int(i_ss_info['ssid']) for i_ss_info in self.wt_ss]
        final_value_ss_id = [int(i_ss_info['ssid']) for i_ss_info in self.ss_info]
        perturbation_names = [i_ss_info['id'] for i_ss_info in self.ss_info]
        dict_fields = ['parameter_name', 'parameter_change', 'parameter_change_percentage', 'parameter_value',
                       'initial_ss', 'final_ss', 'experiment_id', 'acetate']
        experiment_info = dict(zip(dict_fields,
                                   [parameter_name, parameter_change, parameter_change_percentage,
                                    parameter_value, initial_value_ss_id, final_value_ss_id, perturbation_names,
                                    essential_parameter_value]))

        # convert final_ss to dictionary suitable for data frame creation
        concentration_name, concentration_value, sample_name_info = self.create_ss_dict([i_ss["y"]
                                                                                         for i_ss in self.ss_info],
                                                                                        variable_type='metabolite',
                                                                                        noise=self.noise)
        flux_name, flux_value, _ = self.create_ss_dict([i_ss["flux"] for i_ss in self.ss_info], variable_type='flux',
                                                       noise=self.noise)

        # convert experiment_info to dict consistent with concentration_value and flux_value
        for i_field_name in experiment_info:
            new_field_value = self.create_other_value_dict(experiment_info[i_field_name],
                                                           number_of_samples=self.sample_size, noise=self.noise)
            experiment_info[i_field_name] = new_field_value

        experiment_info.update(zip(concentration_name, concentration_value))
        experiment_info.update(zip(flux_name, flux_value))
        experiment_info.update({"sample_name": sample_name_info})

        # prepare list of column names for dataframe
        dict_fields = experiment_info.keys()
        # experiment_info_df = pd.DataFrame(experiment_info, columns=dict_fields)

        import pdb;pdb.set_trace()
        return None

    @staticmethod
    def create_ss_dict(ss_info, variable_type, noise=0):
        """create dictionary of all ss values variables-wise for use in creating data frames"""
        if noise:
            number_variables = len(ss_info[0][0])
        else:
            number_variables = len(ss_info[0])
        variable_name_info = [variable_name(variable_type, j_variable) for j_variable in range(0, number_variables)]

        variable_value_info = []
        if noise:
            for j_variable in range(0, number_variables):
                j_variable_info = []
                for i_sample_id, i_sample_info in enumerate(ss_info):
                    for i_experiment_info in i_sample_info:
                        j_variable_info.append(i_experiment_info[j_variable])
                variable_value_info.append(j_variable_info)

            sample_name_info = []
            for i_sample_id, i_sample_info in enumerate(ss_info):
                for _ in i_sample_info:
                    sample_name_info.append('sample_{}'.format(i_sample_id))
        else:
            for j_variable in range(0, number_variables):
                j_variable_info = []
                for i_experiment_id, i_experiment_info in enumerate(ss_info):
                    j_variable_info.append(i_experiment_info[j_variable])
                variable_value_info.append(j_variable_info)

            sample_name_info = []
            i_sample_id = 0
            for _ in ss_info:
                sample_name_info.append('sample_{}'.format(i_sample_id))

        return variable_name_info, variable_value_info, sample_name_info

    @staticmethod
    def create_other_value_dict(other_info, number_of_samples, noise=0):
        """create dictionary of other values based on number of samples
        to create consistent dict for data frame creation"""
        all_sample_final_ss = []
        for _ in range(0, number_of_samples):
            for i_experiment_value in other_info:
                all_sample_final_ss.append(i_experiment_value)
        return all_sample_final_ss


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
    wt_ss, wt_dynamics = model_1.run_initial_sim([default_parameters], ['default_parameters'])

    # all parameter perturbations
    parameter_perturbation = [{"wt": 0}, {"ac": 1}, {"ac": 4}, {"ac": 9}, {"ac": -.1}, {"ac": -.5},
                              {"k1cat": .1}, {"k1cat": .5}, {"k1cat": 1}, {"k1cat": -.1}, {"k1cat": -.5},
                              {"V3max": .1}, {"V3max": .5}, {"V3max": 1}, {"V3max": -.1}, {"V3max": -.5},
                              {"V2max": .1}, {"V2max": .5}, {"V2max": 1}, {"V2max": -.1}, {"V2max": -.5}]

    experiment_id = ['experiment_{}'.format(parameter_id) for parameter_id, _ in enumerate(parameter_perturbation)]
    experiment_details = model_1.change_parameter_values(parameter_perturbation)

    # call model.simulate to get initial (WT) steady state for all parameter sets strating from same y0
    model_1.sim_model(parameter=experiment_details, experiment_ids=experiment_id, initial_value=wt_ss['y'])
    import pdb;pdb.set_trace()
    print('Done')

    # call model.simulate to get perturbed steady state for all
    # parameter perturbation from corresponding steady states y0
    # sim_perturbation_result = model_1.sim_model

    # call model.change_parameter
