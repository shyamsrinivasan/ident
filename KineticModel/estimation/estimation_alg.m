% kinetic modeling with noisy metabolomics data 
% test using kotte model 
% algorithm :
% min |vmodel - vexpt|
%  st Svmodel = 0;
%     vmodel = f(x,p);
%     p >= pmin;
%     p <= pmax;
% given : x(expt), vexpt

% experimental data is generated through addition of noise to actual kotte model
% vmodel will use convenience kinetics for estimating fluxes

% all parameters available
plist = {'K1ac','K3fdp','L3fdp','K3pep','K2pep','vemax','KeFDP','ne',...
        'd','V4max','k1cat','V3max','V2max','K1pep','K2fdp','rhoA'};

% generate experimental data - get initial ss
tspan = 0:0.1:300;
[xdyn,fdyn,xss1,fss1,opts] = run_nonoise(tspan);

% backup parameters and initial conditions
ival_bkp = opts.x0;
odep_bkp = opts.odep;

% calls to optimize fluxes serially
% use optimized parameters from previous iterations
%% perturbation to flux 1
opts.x0 = xss1;
opts.tspan = 0:.1:300;
opts.odep = odep_bkp;
exp_pid = 11; % 'k1cat'
opts.odep(exp_pid) = 2; % or 2
odep_bkp_3 = opts.odep;
[~,~,xss2_3,fss2_3] = perturb_nonoise(opts);

%% optimize uptake flux parameters - 'K1ac','k1cat','K1pep'
fprintf('\nOptimizing parameters for flux 1.....\n');
[x_opt_3,opt_id_3,~,fval_3] = flux1(opts,xss2_3,fss2_3,plist);

% check if new parameters give the same perturbed flux value
opts.x0 = xss1;
opts.odep = odep_bkp_3;
opts.tspan = 0:.1:50000;
opts.odep(opt_id_3) = x_opt_3;
[~,~,xss4_3,fss4_3] = check_conkin_kotte(opts);

%% perturbation to flux 2
% perutrb model from ss
opts.x0 = xss1;
opts.tspan = 0:.1:300;
exp_pid = 13; % 'V2max'
opts.odep(exp_pid) = .5; 
odep_bkp_1 = opts.odep;
[~,~,xss2_1,fss2_1] = perturb_nonoise(opts);

%% optimize flux2 parameters - 'K2pep','V2max','K2fdp'
fprintf('\nOptimizing parameters for flux 2.....\n');
[x_opt_1,opt_id_1,~,fval_1] = flux2(opts,xss2_1,fss2_1,plist);

% check if new parameters give the same perturbed flux value
opts.x0 = xss1;
opts.odep = odep_bkp_1;
opts.odep(opt_id_1) = x_opt_1;
opts.tspan = 0:.1:500;
[~,~,xss4_1,fss4_1] = check_conkin_kotte(opts);

%% perturbation to flux 3
% restore from backup
opts.x0 = xss1;
opts.tspan = 0:.1:300;
opts.odep = odep_bkp;
exp_pid = 12; % 'V3max'
opts.odep(exp_pid) = 2; % or 2
odep_bkp_2 = opts.odep;
[~,~,xss2_2,fss2_2] = perturb_nonoise(opts);

%% optimize flux 3 with - 'K3fdp','K3pep','V3max','rhoA'
fprintf('\nOptimizing parameters for flux 3.....\n');
[x_opt_2,opt_id_2,~,fval_2] = flux3(opts,xss2_2,fss2_2,plist);

% check if new parameters give the same perturbed flux value
opts.x0 = xss2_2;
opts.odep = odep_bkp_2;
opts.tspan = 0:.1:500;
opts.odep(opt_id_2) = x_opt_2;
[~,~,xss4_2,fss4_2] = check_conkin_kotte(opts);

%% perturbation to flux 4
opts.x0 = xss1;
opts.tspan = 0:.1:300;
opts.odep = odep_bkp;
exp_pid = 10; % 'V4max'
opts.odep(exp_pid) = .8; % or 2
odep_bkp_4 = opts.odep;
[~,~,xss2_4,fss2_4] = perturb_nonoise(opts);

%% optimize parameters for flux 4 - 'V4max'
fprintf('\nOptimizing parameters for flux 4.....\n');
[x_opt_4,opt_id_4,~,fval_4] = flux4(opts,xss2_4,fss2_4,plist);

% check if new parameters give the same perturbed flux value
opts.x0 = xss2_4;
opts.odep = odep_bkp_4;
opts.odep(opt_id_4) = x_opt_4;
opts.tspan = 0:.1:20000;
[~,~,xss4_4,fss4_4] = check_conkin_kotte(opts);

%% perturbation to transcriptional regulation flux
opts.x0 = xss1;
opts.tspan = 0:.1:300;
opts.odep = odep_bkp;
exp_pid = 6; % 'V3max'
opts.odep(exp_pid) = 2; % or 2
odep_bkp_5 = opts.odep;
[~,~,xss2_5,fss2_5] = perturb_nonoise(opts);

%% optimize transcriptional regulation flux
fprintf('\nOptimizing parameters for flux 4.....\n');
[x_opt_5,opt_id_5,~,fval_5] = flux5(opts,xss2_5,fss2_5,plist);

% check if new parameters give the same perturbed flux value
opts.x0 = xss1;
opts.odep = odep_bkp_5;
opts.odep(opt_id_5) = x_opt_5;
opts.tspan = 0:.1:20000;
[~,~,xss4_5,fss4_5] = check_conkin_kotte(opts);

%% check dynamics with new optimized parameters
opts.x0 = xss1;
opts.x0 = ival_bkp;
opts.odep = odep_bkp;
opts.odep(opt_id_1) = x_opt_1;
opts.odep(opt_id_2) = x_opt_2;
opts.odep(opt_id_3) = x_opt_3;
opts.odep(opt_id_4) = x_opt_4;
opts.odep(opt_id_5) = x_opt_5;
opts.tspan = 0:.1:500;
[~,~,xss4_final,fss4_final] = check_conkin_kotte(opts);

% check dynamics with new optimized parameters
% opts.x0 = xss1;
% opts.x0 = ival_bkp;
% opts.odep = odep_bkp;
% opts.odep(opt_id_1) = x_opt_1;
% opts.tspan = 0:.1:300;
% [~,~,xss4_1,fss4_1] = check_conkin_kotte(opts);


% check dynamics with new optimized parameters
% opts.x0 = xss1;
% opts.x0 = ival_bkp;
% opts.odep = odep_bkp;
% opts.odep(opt_id_1) = x_opt_1;
% opts.odep(opt_id_2) = x_opt_2;
% opts.tspan = 0:.1:300;
% [~,~,xss4_2,fss4_2] = check_conkin_kotte(opts);





%% test - generate model data from convinience kinetics
% tspan = 0:0.1:2000;
% [x3dyn,f3dyn,xss3,fss3] = ckintest_script(tspan);
%% run flux parameter estimation
% optimization problem based on minimization of norm of difference in
% fluxes to estimate parameters for one flux
% get expression for objective function and grad(obj) for fluxi from casadi 
% [FX,DFX,D2FX,DFp,fx_sym,x_sym,p_sym,FX_flux] = obj_CAS();



%% noisy stochastic model - uses ode23 to solve
% tspan = 0:0.1:500;
% run_withnoise

% perturb noisy model from ss
% perturb enzyme 2 (flux(4),vEX)
% odep.p(14) = .5; % 2;
% perturb_noise


% optimize for parameters by minimizing L2 norm of fluxes






