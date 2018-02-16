% run ode to simulate data
solver_options.tspan = 0:.1:250;
solver_options.x0 = [5;1;1];
solver_options.reltol = 1e-6;
solver_options.abstol = 1e-8;
fun_p = [.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, .1];
solution = run_ode(@kotteode, solver_options, fun_p, @kotteflux);
% plot results for verification
subplot(211)
plot(solution.t, solution.y)
subplot(212)
plot(solution.t, solution.flux)

% parameter perturbations
perturbations = perturbations_info();
solver_options.x0 = solution.yss;
perturbed_solution = perturb_parameters(@kotteode, solver_options,...
                                        fun_p, perturbations, @kotteflux);
                                    
% arrange simulated data for parameter estimation
 
% use simulated data from 3 different perturbations to solve nlae for
% parameters



