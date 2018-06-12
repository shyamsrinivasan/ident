from create_experiment_data import retrieve_experimental_data_from_file
import pandas as pd
import os.path


# get identifiability info all fluxes from ident files
v1_file_name = os.path.join(os.getcwd(), 'ident/ident_v1_kcat')
ident_index_label = ['sample_name', 'data_set_id']
# retrieve identifiability info from file
ident_df = retrieve_experimental_data_from_file(v1_file_name, ident_index_label)
# all_parameter_info = process_ident(ident_df, arranged_data_df)
idx = pd.IndexSlice

all_parameter_names = ident_df["parameter_name"].unique().tolist()
number_experiments = len(all_parameter_names)
exp_column_ids = ['experiment_{}_id'.format(i_experiment) for i_experiment in range(0, number_experiments)]
for i_parameter, i_parameter_name in enumerate(all_parameter_names):
    # get all data sets identifying each parameter
    identifying_df = ident_df[(ident_df["parameter_name"] == i_parameter_name) & (ident_df["identified"])]
    for i_experiment, i_experiment_pos in enumerate(exp_column_ids):
        exp_frequency = identifying_df[i_experiment_pos].value_counts()
        pass
    pass

print('Run Complete\n')