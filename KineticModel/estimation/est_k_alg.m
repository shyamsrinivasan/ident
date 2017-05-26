% estimate parameters for same model (kotte) after perturbations
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
[x_opt_3,opt_id_3,~,fval_3] = flux1_k(opts,xss2_3,fss2_3,plist);

%% check if new parameters give the same perturbed flux value
opts.x0 = xss1;
opts.odep = odep_bkp_3;
opts.tspan = 0:.1:500;
opts.odep(opt_id_3) = x_opt_3;
[~,~,xss4_3,fss4_3] = check_kin_kotte(opts);

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
[x_opt_1,opt_id_1,~,fval_1] = flux2_k(opts,xss2_1,fss2_1,plist);

%% check if new parameters give the same perturbed flux value
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
[x_opt_2,opt_id_2,~,fval_2] = flux3_k(opts,xss2_2,fss2_2,plist);

%% check if new parameters give the same perturbed flux value
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
[x_opt_4,opt_id_4,~,fval_4] = flux4_k(opts,xss2_4,fss2_4,plist);

%% check if new parameters give the same perturbed flux value
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

%% check if new parameters give the same perturbed flux value
opts.x0 = xss1;
opts.odep = odep_bkp_5;
opts.odep(opt_id_5) = x_opt_5;
opts.tspan = 0:.1:20000;
[~,~,xss4_5,fss4_5] = check_conkin_kotte(opts);
