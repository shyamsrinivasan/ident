% run ode to simulate data
solver_options.tspan = 0:.1:250;
solver_options.x0 = [10;5;1];
solver_options.reltol = 1e-6;
solver_options.abstol = 1e-8;
fun_p = [.1, .1, 4e6, .1, .3, 1.1, .45, 2, .25, .2, 1, 1, 1, .1];
solution = run_ode(@kotteCKode, solver_opts, fun_p);
% use simulated data from 3 different perturbations to solve nlae for
% parameters