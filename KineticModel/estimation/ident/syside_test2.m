% generate experimental data - get initial ss
tspan = 0:0.1:100;
[xdyn,fdyn,xss,fss,opts] = run_nonoise(tspan);

% backup parameters and initial conditions
ival_bkp = opts.x0;
odep_bkp = opts.odep;

% perturb system from non-noisy initial conditions
pt_sol_id = [1 2 3];
[exp_sol,no_noise_sol] = dopert_nonoise_dyn(opts,xss,fss,odep_bkp,pt_sol_id);
close all

% new sysid_test with functions
[xdyn_fun,ode,x,p] = RK4integrator_cas(@kotte_CAS);

xinit = no_noise_sol(4).xss;
xexp_data = no_noise_sol(1).xdyn(:,2:end);

x_sym = xdyn_fun(no_noise_sol(4).xss,repmat(p,1,3000));
x_error = (xexp_data-x_sym);
obj = .5*dot(x_error,x_error);

% setup nlp
nlp = struct('x',p,'f',obj);
solver = casadi.nlpsol('solver','ipopt',nlp);
sol_opt = solver('x0',opts.odep(2:end)');