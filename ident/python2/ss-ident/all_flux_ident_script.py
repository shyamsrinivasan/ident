from create_experiment_data import retrieve_experimental_data_from_file
import pandas as pd
import os.path


# get identifiability info all fluxes from ident files
v1_file_name = os.path.join(os.getcwd(), 'ident/ident_v1_kcat')
ident_index_label = ['sample_name', 'data_set_id']
# retrieve identifiability info from file
ident_df = retrieve_experimental_data_from_file(v1_file_name, ident_index_label)
# original set of experiments
original_experiment_file = os.path.join(os.getcwd(), 'exp/experiments')
exp_df = retrieve_experimental_data_from_file(data_file_name=original_experiment_file,
                                              multi_index_label=['sample_name', 'experiment_id'])
all_possible_experiments = exp_df.index.levels[1].tolist()
# all_parameter_info = process_ident(ident_df, arranged_data_df)
idx = pd.IndexSlice

all_parameter_names = ident_df["parameter_name"].unique().tolist()
number_experiments = len(all_parameter_names)
exp_column_ids = ['experiment_{}_id'.format(i_experiment) for i_experiment in range(0, number_experiments)]
for i_parameter, i_parameter_name in enumerate(all_parameter_names):
    # get all data sets identifying each parameter
    all_experiment_pos_info = []
    identifying_df = ident_df[(ident_df["parameter_name"] == i_parameter_name) & (ident_df["identified"])]
    for i_experiment, i_experiment_pos in enumerate(exp_column_ids):
        exp_frequency = identifying_df[i_experiment_pos].value_counts()
        available_experiments = exp_frequency.index.tolist()
        missing_experiments = [i_exp for i_exp in all_possible_experiments if i_exp not in available_experiments]
        missing_experiment_value = [0] * len(missing_experiments)
        missing_series = pd.Series(missing_experiment_value, index=missing_experiments)
        new_experiment_frequency = exp_frequency.append(missing_series).rename(exp_frequency.name)
        all_experiment_pos_info.append(new_experiment_frequency)
        pass
    all_keys = [j_series.name for j_series in all_experiment_pos_info]
    new_df = pd.concat(all_experiment_pos_info,
                       keys=all_keys,
                       names=['experiment_pos', 'experiment_id']).to_frame(name='frequency').reset_index()
    complete_df = new_df.assign(parameter_name=[i_parameter_name] * len(new_df.index),
                                flux_name=identifying_df["flux_name"].unique().tolist() * len(new_df.index))

    pass

print('Run Complete\n')
