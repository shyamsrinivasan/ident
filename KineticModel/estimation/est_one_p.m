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

%% optimize flux2 parameters - 'K2pep','V2max','K2fdp'
fprintf('\nOptimizing parameters for flux 2.....\n');
[x_opt_1,opt_id_1,~,fval_1] = flux2_k(opts,xss2_3,fss2_3,plist);

%% optimize flux 3 with - 'K3fdp','K3pep','V3max','rhoA'
fprintf('\nOptimizing parameters for flux 3.....\n');
% [x_opt_2,opt_id_2,~,fval_2] = flux3_k(opts,xss2_2,fss2_2,plist);
% testing the use of scip as a global optimizer for nlp
[x_opt_2,opt_id_2,~,fval_2] = flux3_k(opts,xss2_3,fss2_3,plist);

%% optimize parameters for flux 4 - 'V4max'
fprintf('\nOptimizing parameters for flux 4.....\n');
[x_opt_4,opt_id_4,~,fval_4] = flux4_k(opts,xss2_3,fss2_3,plist);
