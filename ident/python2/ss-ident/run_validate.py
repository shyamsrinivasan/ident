from run_sims import ModelSim


class ValidateSim(ModelSim):
    def __init__(self, rhs_fun, flux_fun, noise=0, **kwargs):
        super(ValidateSim, self).__init__(rhs_fun, flux_fun, noise, **kwargs)
        # self.rhs_fun = rhs_fun
        # self.flux_fun = flux_fun
        # self.noise = noise

        # all parameter perturbations to be tested
        try:
            self.test_perturbations = kwargs['test_perturbations']
        except KeyError:
            self.test_perturbations = [{"wt": 0}, {"ac": 1}, {"ac": 4}, {"ac": 9}, {"ac": -.1}, {"ac": -.5},
                                       {"k1cat": .1}, {"k1cat": .5}, {"k1cat": 1}, {"k1cat": -.1}, {"k1cat": -.5},
                                       {"V3max": .1}, {"V3max": .5}, {"V3max": 1}, {"V3max": -.1}, {"V3max": -.5},
                                       {"V2max": .1}, {"V2max": .5}, {"V2max": 1}, {"V2max": -.1}, {"V2max": -.5}]

        # self.wt_ss = []
        # self.wt_dynamic = []
        # self.perturbation_ss = []
        # self.perturbation_dynamic = []
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
        # data_set_id = []
        estimate_data_set_info = []
        for i_index, i_value in enumerate(parameter_list):
            # data_set_id.append(estimate_info['data_sets'][i_index])
            # estimate_id.append('estimate_{}'.format(i_index))
            estimate_data_set_info.append(('estimate_{}'.format(i_index), estimate_info['data_sets'][i_index][0],
                                           estimate_info['data_sets'][i_index][1]))
            for i_key in parameter_name_value_pair[i_index].keys():
                i_value[i_key] = parameter_name_value_pair[i_index][i_key]

        return parameter_list, estimate_data_set_info

    @staticmethod
    def separate_initial_perturbation(all_sim_results):
        """separate initial ss simulations from perturbation simulation
        for further processing of perturbation sims"""

        initial_sims = []
        perturbation_sims = []

        for j_result in all_sim_results:
            if j_result['initial']:
                initial_sims.append(j_result)
            elif j_result['perturbation']:
                perturbation_sims.append(j_result)

        return initial_sims, perturbation_sims

    @staticmethod
    def separate_ss_dyn(all_results):
        """get ss and dynamic information separately"""

        ss_info = []
        dynamic_info = []

        for j_result in all_results:
            ss_info.append(j_result['ss'])
            dynamic_info.append(j_result['dynamic'])

        return ss_info, dynamic_info

    def create_ss_perturbation_dict(self, all_results):

        # separate initial simulation values from perturbation simulation values
        initial_sims, perturbation_sims = self.separate_initial_perturbation(all_results)

        # get initial ss values only
        initial_ss, _ = self.separate_ss_dyn(initial_sims)

        # get perturbation ss values only
        perturbation_ss, _ = self.separate_ss_dyn(perturbation_sims)

        # create dict of perturbation_ss values for writing to df and file
        import pdb;pdb.set_trace()
        print('What to do next?\n')






    # def run_initial_sim(self, parameter, parameter_ids=(), **kwargs):
    #     import pdb;pdb.set_trace()
    #     wt_ss, wt_dynamics = super(ValidateSim, self).run_initial_sim(parameter, parameter_ids, **kwargs)
    #     return wt_ss, wt_dynamics
    #
    # def validate_model(self, parameter_estimates, estimate_info):
    #     """run parallel validation method"""
    #     # run parallel initial sim (based on setup parallel ode for multiple parameter sets)
    #     import pdb;pdb.set_trace()
    #     estimate_ids = [j_value[0] for j_value in estimate_info]
    #     initial_ss, initial_dynamics = self.run_initial_sim(parameter_estimates, estimate_ids)
    #
    #     import pdb;pdb.set_trace()
    #     sim_result = setup_parallel_ode(ode_rhs_fun=self.rhs_fun, flux_fun=self.flux_fun,
    #                                     parameters=parameter_estimates, y0=self.wt_y0, t_final=self.t_final,
    #                                     experiment_id=estimate_ids, ode_opts=self.ode_opts, i_value_opt=0,
    #                                     parameter_opt=1)
    #     # run parallel perturbation sim for each parameter estimate
    #     import pdb;pdb.set_trace()
    #     return sim_result
