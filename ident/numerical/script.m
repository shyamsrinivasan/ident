% run ode to simulate data
solver_options.tspan = 0:.1:250;
solver_options.x0 = [5;1;1];
solver_options.reltol = 1e-6;
solver_options.abstol = 1e-8;
fun_p = [.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, .1];
solution = run_ode(@kotteode, solver_options, fun_p, @kotteflux);
% plot results for verification
plot2d_dynamic(solution, 1, 1);

% parameter perturbations
perturbations = perturbations_info();
solver_options.x0 = solution.yss;
perturbed_solution = perturb_parameters(@kotteode, solver_options,...
                                        fun_p, perturbations, @kotteflux);
plot2d_dynamic(perturbed_solution, 1, 1);
                                    
% arrange simulated data for parameter estimation
arrange_experimental_data(perturbed_solution, 3, 1:5);
 
% use simulated data from 3 different perturbations to solve nlae for
% parameters



