from names_strings import true_parameter_values
from run_ident import ModelIdent
from run_validation import ValidateSim
from mpi4py import MPI
import numpy as np
import kotte_model
import os.path


def main():
    name = MPI.Get_processor_name()
    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()
    print('I am  %s rank %d (total %d)' % (name, rank, size))

    if rank == 0:
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

        import pdb;pdb.set_trace()
        v1_ident.process_ident()

        v1_ident.get_parameter_value()

        # validate estimated parameter values
        user_ode_opts = {'iter': 'Newton', 'discr': 'Adams', 'atol': 1e-10, 'rtol': 1e-10,
                         'time_points': 200, 'display_progress': True, 'verbosity': 30}
        # initial ss to begin all simulations from
        y0 = np.array([5, 1, 1])
        # get and set true parameter values, if available separately
        default_parameters = true_parameter_values()

        validate_obj = ValidateSim(kotte_model.kotte_ck_ode, kotte_model.kotte_ck_flux, **{'kinetics': 2,
                                                                                           'ode_opts': user_ode_opts,
                                                                                           't_final': 200,
                                                                                           'wt_y0': y0,
                                                                                           'i_parameter': default_parameters,
                                                                                           'sample_size': 1,
                                                                                           'noise_std': 0.05})

        parameter_estimates, estimate_info = validate_obj.create_parameter_list(v1_ident.select_values)

        import pdb;pdb.set_trace()
        validate_obj.validate_model(parameter_estimates, estimate_info=estimate_info)

        # get parameter value plot
        v1_ident.parameter_values_plot(default_parameter_values, violin=True, box=False, bins=1)

        v1_ident.identifiability_plot()

        v1_ident.exp_info_plot()

        import pdb;pdb.set_trace()
        print('Done\n')


if __name__ == '__main__':
    main()
