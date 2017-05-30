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
exp_pval = [.2;.5;.8;1.2;1.5;1.8;2];
[xss_1,fss_1] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

%% perturbation to flux 2
% perutrb model from ss
opts.x0 = xss1;
opts.tspan = 0:.1:300;
opts.odep = odep_bkp;
exp_pid = 13; % 'V2max'
exp_pval = [.2;.5;.8;1.2;1.5;1.8;2];
[xss_2,fss_2] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

%% perturbation to flux 3
% restore from backup
opts.x0 = xss1;
opts.tspan = 0:.1:300;
opts.odep = odep_bkp;
exp_pid = 12; % 'V3max'
exp_pval = [.2;.5;.8;1.2;1.5;1.8;2];
[xss_3,fss_3] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

%% perturbation to flux 4
opts.x0 = xss1;
opts.tspan = 0:.1:10000;
opts.odep = odep_bkp;
exp_pid = 10; % 'V4max'
exp_val = [0;.1;.4;.5;.6;.8;1];
[xss_4,fss_4] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

%% optimize flux 3 with - 'K3fdp','K3pep','V3max','rhoA'
opts.odep = odep_bkp;
fprintf('\nOptimizing parameters for flux 3.....\n');
[x_opt_1,opt_id_1,~,fval_1] =...
flux3_k(opts,[xss_1,xss_2,xss_3,xss_4],[fss_1,fss_2,fss_3,fss_4],plist);

fprintf('\nOptimizing parameters for flux 3.....\n');
[x_opt_2,opt_id_2,~,fval_2] =...
flux3_k(opts,[xss_1,xss_2,xss_3],[fss_1,fss_2,fss_3],plist);

fprintf('\nOptimizing parameters for flux 3.....\n');
[x_opt_3,opt_id_3,~,fval_3] =...
flux3_k(opts,[xss_1,xss_2],[fss_1,fss_2],plist);

fprintf('\nOptimizing parameters for flux 3.....\n');
[x_opt_4,opt_id_4,~,fval_4] =...
flux3_k(opts,xss_1,fss_1,plist);

%% check if new parameters give the same perturbed flux value
% rerun perturbations with new parameters
opts.x0 = xss1;
opts.odep = odep_bkp;
opts.tspan = 0:.1:10000;
opts.odep(opt_id_1) = x_opt_1;
exp_pid = 11;
exp_pval = [.2;.5;.8;1.2;1.5;1.8;2];
[xss4_1,fss4_1] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);

% [~,~,xss4_1,fss4_1] = check_kin_kotte(opts);



