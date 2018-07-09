from names_strings import true_parameter_values
from run_ident import ModelIdent
import kotte_model
import os.path

# extract experimental data from file
default_parameter_values = true_parameter_values()

v1_ident = ModelIdent(ident_fun=kotte_model.flux_1_kcat_ident,
                      arranged_data_file_name=os.path.join(os.getcwd(), 'exp/exp_v1_2_experiments'),
                      ident_data_file_name=os.path.join(os.getcwd(), 'ident/ident_v1_kcat'),
                      **{'original_exp_file': os.path.join(os.getcwd(), 'exp/experiments'),
                         'flux_id': 1, 'flux_choice': 2,
                         'values_figure': os.path.join(os.getcwd(), 'results/v1_kcat_parameter_values.eps'),
                         'ident_figure': os.path.join(os.getcwd(), 'results/v1_kcat_ident.eps'),
                         'exp_figure': os.path.join(os.getcwd(), 'results/v1_kcat_exp.eps'),
                         'figure_format': 'eps'})
# test identifiability
print('Practical Identifiability Analysis of v1 with 2 parameters: k1cat and K1ac\n')
ident_data_df = v1_ident.perform_ident()

v1_ident.process_ident()

v1_ident.get_parameter_value()

# validate estimated parameter values


# get parameter value plot
v1_ident.parameter_values_plot(default_parameter_values, violin=True, box=False, bins=1)

v1_ident.identifiability_plot()

v1_ident.exp_info_plot()

import pdb;pdb.set_trace()
print('Done\n')
