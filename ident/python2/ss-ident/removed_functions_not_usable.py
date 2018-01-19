def experiments_in_ident_data(data_p_boolean, data_exp, exp_types, data_id):
    """get all data sets that identify parameter j and get different types of experiments in these data sets"""
    number_of_experiments_per_data = 3
    number_exp_types = len(exp_types)
    data_identifying_p = []
    exp_data_parameter_info = []
    for j_p, data_set_boolean in enumerate(data_p_boolean):
        data_identifying_parameter_j = [j_data for j_data, value in enumerate(data_set_boolean) if value]
        data_identifying_p.append(data_identifying_parameter_j)
        chosen_data_experiments = [j_data_set for data_set_id, j_data_set in enumerate(data_exp["experiment_id"])
                                   if data_set_id in data_identifying_parameter_j]
        chosen_data_experiments = np.array(chosen_data_experiments)
        # get number of experiments of each type at each position
        exp_types_classification = []
        for i_exp_type in exp_types:
            try:
                pos_based_lst = [[chosen_data_experiments[:, iexp] == np.array(j_exp_id) for j_exp_id in i_exp_type]
                                 for iexp in range(0, number_of_experiments_per_data)]
            except IndexError:
                pos_based_lst = []
            exp_types_classification.append(pos_based_lst)
        # get experiment type numbers and percentages within data set at each position
        total_identifying_experiments = chosen_data_experiments.shape[0]
        counter_exp_type = 0
        # print('Parameter {}'.format(j_p))
        exp_type_all_pos_info = []
        for j_exp_type, exp_id in zip(exp_types_classification, exp_types):
            all_exp_types_info = []
            for i_pos, i_exp_pos in enumerate(j_exp_type):
                all_pos_info = []
                # print('Experiment Position in Data sets {}'.format(i_pos))
                for k_exp_at_i_exp_pos_in_j_exp_type, k_exp in zip(i_exp_pos, exp_id):
                    total_number_of_k_exp_at_i_pos_in_j_exp_type = sum(k_exp_at_i_exp_pos_in_j_exp_type)
                    percentage_of_k_exp_at_i_pos_in_j_exp_type = \
                        float(sum(k_exp_at_i_exp_pos_in_j_exp_type))/float(total_identifying_experiments)
                    info = {'type':counter_exp_type, 'id':k_exp,
                            'occurrence':total_number_of_k_exp_at_i_pos_in_j_exp_type,
                            'occurrence percentage':percentage_of_k_exp_at_i_pos_in_j_exp_type,
                            'position':i_pos}
                    all_pos_info.append(info)
                    pass
                # all_pos_info.append(all_pos_info)
                all_exp_types_info.append(all_pos_info)
            counter_exp_type += 1
            exp_type_all_pos_info.append(all_exp_types_info)
        exp_data_parameter_info.append(exp_type_all_pos_info)
    return exp_data_parameter_info