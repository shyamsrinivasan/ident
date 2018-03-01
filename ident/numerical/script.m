% specify ode fun and flux fun
odefun = @kotteode;
fluxfun = @kotteflux;

% run ode to simulate data
solver_options.tspan = 0:.1:250;
solver_options.x0 = [5;1;1];
solver_options.reltol = 1e-6;
solver_options.abstol = 1e-8;
fun_p = [.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, .1];
solution = run_ode(odefun, solver_options, fun_p, fluxfun);
% plot results for verification
plot2d_dynamic(solution, 1, 1);

% parameter perturbations
perturbations = perturbations_info();
solver_options.x0 = solution.yss;
perturbed_solution = perturb_parameters(odefun, solver_options,...
                                        fun_p, perturbations, fluxfun);
plot2d_dynamic(perturbed_solution, 1, 1);
                                    
% arrange simulated data for parameter estimation
nexpts = 3;
flux_id = [1, 4, 3, 5];
% choose only experiments w/0 v3max perturbation
selected_expts = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 17, 18, 19, 20, 21];
experimental_data =...
    arrange_experimental_data(perturbed_solution,...
                              nexpts,...
                              flux_id,...
                              selected_expts);

% run numerical parameter estimation (w/MWC model) on all combinations 
% in experimental_data 
ident.fun = @estimate_v3_parameter;
ident.initial_value = [1;.1;.1];
ident.typical_value = [1;1e-1;1e-1];
solution = do_numerical_ident(experimental_data,ident);

% use simulated data from 3 different perturbations to solve nlae for
% parameters
% initial_p_val = [1;.1;.1];
% typical_p_val = [1;1e-1;1e-1];
% nlae_solution(@estimate_v3_CK_parameter, experimental_data(1, :), initial_p_val, typical_p_val);

% estimate parameters for MWC model
% initial_p_val = [1;.1;.1];
% typical_p_val = [1;1e-1;1e-1];
% nlae_solution(@estimate_v3_parameter, experimental_data(1, :), initial_p_val, typical_p_val);





