% generate experimental data - get initial ss
tspan = 0:.1:300;
[xdyn,fdyn,xss,fss,opts] = run_nonoise(tspan);

% backup parameters and initial conditions
ival_bkp = opts.x0;
odep_bkp = opts.odep;

% perturb system from non-noisy initial conditions
opts.tspan = 0:.1:200;
pt_sol_id = [1 2 3];
[exp_sol,no_noise_sol] = dopert_nonoise_dyn(opts,xss,fss,odep_bkp,pt_sol_id);
close all

optim_opts = struct('pname','K1ac','nc',3,'nf',6,...
                    'odep',odep_bkp,...
                    'tspan',opts.tspan,...
                    'x0',xss);

% new sysid_test with functions
[xdyn_fun,ode,x,p,ident_c,p_useless,acetate] =...
RK4integrator_cas(@kotteCASident,optim_opts);

xinit = no_noise_sol(4).xss;
xexp_data = no_noise_sol(1).xdyn(:,1:100:2001);

npts = length(opts.tspan)-1;
x_sym =...
xdyn_fun(no_noise_sol(4).xss,repmat(p,1,npts),.1,repmat(opts.odep(14:16)',1,npts),.1);
% add initial value
x_sym = [casadi.DM(no_noise_sol(4).xss) x_sym];
x_model_sym = x_sym(:,1:100:2001);
x_error = (xexp_data-x_model_sym);
obj = .5*dot(x_error,x_error);

% setup nlp
nlp = struct('x',p,'f',obj);
solver = casadi.nlpsol('solver','ipopt',nlp);
p0 = opts.odep(2:13)';
p0(6) = .1;
sol_opt = solver('x0',opts.odep(2:13)');

