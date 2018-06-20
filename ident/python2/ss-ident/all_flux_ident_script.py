# from create_experiment_data import retrieve_experimental_data_from_file
from process_exp_details import exp_design_info
from process_exp_details import logical_values
from plot_ident_results import plot_exp_details
import os.path


# get all info on original experiments
original_experiment_file = os.path.join(os.getcwd(), 'exp/experiments')

# get identifiability info all fluxes from ident files
v1_file_name = os.path.join(os.getcwd(), 'ident/ident_v1_kcat')
v2_file_name = os.path.join(os.getcwd(), 'ident/ident_v2')
v3_file_name = os.path.join(os.getcwd(), 'ident/ident_v3_root_1')
v5_file_name = os.path.join(os.getcwd(), 'ident/ident_v5_root_2')
exp_2_file_name_list = [v1_file_name, v2_file_name, v5_file_name]

write_to_file_name = os.path.join(os.getcwd(), 'ident/ident_2_experiments')

df_2 = exp_design_info(list_of_files=exp_2_file_name_list, original_experiment_file=original_experiment_file,
                       write_to_file_name=write_to_file_name, max_number_experiments=2)
plot_exp_details(df_2, max_number_experiments=2, color_bar=True, set_palette=False)

# plot only logical only (experiment is involved/not involved)
logical_df = df_2.applymap(logical_values)
plot_exp_details(logical_df, max_number_experiments=2)


# idx = pd.IndexSlice
# new_df = df.loc[idx[:, :], idx['experiment_0_id', :]]
# iris = sns.load_dataset('iris')

exp_3_file_name_list = [v3_file_name]
write_to_file_name = os.path.join(os.getcwd(), 'ident/ident_3_experiments')
df_3 = exp_design_info(list_of_files=exp_3_file_name_list, original_experiment_file=original_experiment_file,
                       write_to_file_name=write_to_file_name, max_number_experiments=3)
plot_exp_details(df_3, max_number_experiments=3, color_bar=True, set_palette=False)

# plot only logical only (experiment is involved/not involved)
logical_df = df_3.applymap(logical_values)
plot_exp_details(logical_df, max_number_experiments=3)

print('Run Complete\n')
