# from create_experiment_data import retrieve_experimental_data_from_file
from process_exp_details import exp_design_info, logical_values
from plot_ident_results import plot_exp_details
import os.path


# get all info on original experiments
original_experiment_file = os.path.join(os.getcwd(), 'exp/experiments')

# get identifiability info all fluxes from ident files
v1_file_name = os.path.join(os.getcwd(), 'ident/ident_v1_kcat')
v2_file_name = os.path.join(os.getcwd(), 'ident/ident_v2')
v3_file_name = os.path.join(os.getcwd(), 'ident/ident_v3')
v5_file_name = os.path.join(os.getcwd(), 'ident/ident_v5')
file_name_list = [v1_file_name, v2_file_name, v3_file_name, v5_file_name]

write_to_file_name = os.path.join(os.getcwd(), 'ident/ident_experiments')

df = exp_design_info(list_of_files=file_name_list, original_experiment_file=original_experiment_file,
                     write_to_file_name=write_to_file_name, max_number_experiments=3)

# idx = pd.IndexSlice
# new_df = df.loc[idx[:, :], idx['experiment_0_id', :]]
# iris = sns.load_dataset('iris')


logical_df = df.applymap(logical_values)
plot_exp_details(logical_df)
plot_exp_details(df, color_bar=True, set_palette=False)

print('Run Complete\n')
