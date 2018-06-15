# from create_experiment_data import retrieve_experimental_data_from_file
from process_exp_details import exp_design_info
from process_exp_details import get_original_experiments
import os.path


# get all info on original experiments
original_experiment_file = os.path.join(os.getcwd(), 'exp/experiments')

# get identifiability info all fluxes from ident files
v1_file_name = os.path.join(os.getcwd(), 'ident/ident_v1_kcat')
v2_file_name = os.path.join(os.getcwd(), 'ident/ident_v2')
v3_file_name = os.path.join(os.getcwd(), 'ident/ident_v3_root_1')
v5_file_name = os.path.join(os.getcwd(), 'ident/ident_v5_root_2')
file_name_list = [v1_file_name, v2_file_name, v3_file_name, v5_file_name]

exp_design_info(list_of_files=file_name_list, original_experiment_file=original_experiment_file,
                max_number_experiments=3)

print('Run Complete\n')
