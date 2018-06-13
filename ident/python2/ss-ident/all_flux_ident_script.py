from create_experiment_data import retrieve_experimental_data_from_file
from collections import defaultdict
import pandas as pd
import itertools as it
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
empty_dict = {}
all_info_dict = defaultdict(list)

# code for adding data to existing dict
# for key, value in it.chain(empty_dict.items(), new_dict.items()):
#     for i_value in value:
#         all_info_dict[key].append(i_value)

all_parameter_names = ident_df["parameter_name"].unique().tolist()
all_parameter_info = []
all_parameter_flux = []
number_experiments = len(all_parameter_names)
exp_column_ids = ['experiment_{}_id'.format(i_experiment) for i_experiment in range(0, number_experiments)]
for i_parameter, i_parameter_name in enumerate(all_parameter_names):
    # get all data sets identifying each parameter
    all_experiment_info = defaultdict(list)
    identifying_df = ident_df[(ident_df["parameter_name"] == i_parameter_name) & (ident_df["identified"])]
    for i_experiment, i_experiment_pos in enumerate(exp_column_ids):
        exp_frequency = identifying_df[i_experiment_pos].value_counts()
        available_experiments = exp_frequency.index.tolist()
        missing_experiments = [i_exp for i_exp in all_possible_experiments if i_exp not in available_experiments]
        missing_experiment_value = [0] * len(missing_experiments)
        missing_series = pd.Series(missing_experiment_value, index=missing_experiments)
        new_experiment_frequency = exp_frequency.append(missing_series).rename(exp_frequency.name)
        experiment_ids = new_experiment_frequency.index.tolist()
        experiment_frequency = [int(i_value) for i_value in new_experiment_frequency.values.tolist()]
        experiment_pos = [i_experiment_pos] * len(experiment_ids)
        # flux_names = identifying_df["flux_name"].unique().tolist()
        new_dict = dict(zip(['experiment_id', 'frequency', 'experiment_pos_id'],
                            [experiment_ids, experiment_frequency, experiment_pos]))
        for key, value in it.chain(empty_dict.items(), new_dict.items()):
            for i_value in value:
                all_experiment_info[key].append(i_value)
    all_parameter_info.append(all_experiment_info)
    all_parameter_flux.append(identifying_df["flux_name"].unique()[0])
    pass
all_frequency = [i_parameter["frequency"] for i_parameter in all_parameter_info]
col_ind_tuple = zip(all_parameter_info[0]['experiment_pos_id'], all_parameter_info[0]['experiment_id'])
col_ind = pd.MultiIndex.from_tuples(col_ind_tuple, names=['experiment_pos_id', 'experiment_id'])
flux_tuple = zip(all_parameter_flux, all_parameter_names)
row_ind = pd.MultiIndex.from_tuples(flux_tuple, names=['flux_name', 'parameter_name'])
df = pd.DataFrame(all_frequency, index=row_ind, columns=col_ind)




print('Run Complete\n')
