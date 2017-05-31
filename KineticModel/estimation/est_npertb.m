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
exp_pval = [.1;.5;1.0;1.5;2];
[xss_1,fss_1] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

%% perturbation to flux 2
% perutrb model from ss
opts.x0 = xss1;
opts.tspan = 0:.1:300;
opts.odep = odep_bkp;
exp_pid = 13; % 'V2max'
exp_pval = [.1;.5;1.0;1.5;2];
[xss_2,fss_2] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

%% perturbation to flux 3
% restore from backup
opts.x0 = xss1;
opts.tspan = 0:.1:300;
opts.odep = odep_bkp;
exp_pid = 12; % 'V3max'
exp_pval = [.1;.5;1.0;1.5;2];
[xss_3,fss_3] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

%% perturbation to flux 4
opts.x0 = xss1;
opts.tspan = 0:.1:10000;
opts.odep = odep_bkp;
exp_pid = 10; % 'V4max'
exp_pval = [0;.1;.3;.5;.7;.9;1];
[xss_4,fss_4] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

%% use only 1 perturbation set flux to get parameters
x_opt_1 = optimize_p(opts,xss_1,fss_1,plist,odep_bkp);

%% use 2 perturbation sets 
x_opt_2 = optimize_p(opts,[xss_1,xss_2],[fss_1,fss_2],plist,odep_bkp);

%% use 3 perturbation sets
x_opt_3 =...
optimize_p(opts,[xss_1,xss_2,xss_3],[fss_1,fss_2,fss_3],plist,odep_bkp);

%% use 4 perturbation sets
x_opt_4 =...
optimize_p(opts,[xss_1,xss_2,xss_3,xss_4],[fss_1,fss_2,fss_3,fss_4],plist,odep_bkp);

%% rerun all perturbations with new parameters from x_opt_1
% check if new parameters give the same perturbed flux value
opts.x0 = xss1;
opts.odep = odep_bkp;
opts.tspan = 0:.1:10000;
allopt_id_1 = cat(2,x_opt_1(:).opt_id);
allx_opt_1 = cat(1,x_opt_1(:).x_opt);
opts.odep(allopt_id_1) = allx_opt_1;

% perturbation to flux 1
exp_pid = 11; % 'k1cat'
exp_pval = [.1;.5;1.0;1.5;2];
[xss4_1,fss4_1] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

% perturbation to flux 2
exp_pid = 13; % 'V2max'
exp_pval = [.1;.5;1.0;1.5;2];
[xss4_2,fss4_2] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

% perturbation to flux 3
% restore from backup
exp_pid = 12; % 'V3max'
exp_pval = [.1;.5;1.0;1.5;2];
[xss4_3,fss4_3] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

% perturbation to flux 4
exp_pid = 10; % 'V4max'
exp_pval = [0;.1;.3;.5;.7;.9;1];
[xss4_4,fss4_4] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

%% rerun all perturbations with new parameters from x_opt_2
% check if new parameters give the same perturbed flux value
opts.x0 = xss1;
opts.odep = odep_bkp;
opts.tspan = 0:.1:10000;
allopt_id_2 = cat(2,x_opt_2(:).opt_id);
allx_opt_2 = cat(1,x_opt_2(:).x_opt);
opts.odep(allopt_id_2) = allx_opt_2;

% perturbation to flux 1
exp_pid = 11; % 'k1cat'
exp_pval = [.1;.5;1.0;1.5;2];
[xss5_1,fss5_1] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

% perturbation to flux 2
exp_pid = 13; % 'V2max'
exp_pval = [.1;.5;1.0;1.5;2];
[xss5_2,fss5_2] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

% perturbation to flux 3
% restore from backup
exp_pid = 12; % 'V3max'
exp_pval = [.1;.5;1.0;1.5;2];
[xss5_3,fss5_3] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

% perturbation to flux 4
exp_pid = 10; % 'V4max'
exp_pval = [0;.1;.3;.5;.7;.9;1];
[xss5_4,fss5_4] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

%% rerun all perturbations with new parameters from x_opt_3
% check if new parameters give the same perturbed flux value
opts.x0 = xss1;
opts.odep = odep_bkp;
opts.tspan = 0:.1:10000;
allopt_id_3 = cat(2,x_opt_3(:).opt_id);
allx_opt_3 = cat(1,x_opt_3(:).x_opt);
opts.odep(allopt_id_3) = allx_opt_3;

% perturbation to flux 1
exp_pid = 11; % 'k1cat'
exp_pval = [.1;.5;1.0;1.5;2];
[xss6_1,fss6_1] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

% perturbation to flux 2
exp_pid = 13; % 'V2max'
exp_pval = [.1;.5;1.0;1.5;2];
[xss6_2,fss6_2] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

% perturbation to flux 3
% restore from backup
exp_pid = 12; % 'V3max'
exp_pval = [.1;.5;1.0;1.5;2];
[xss6_3,fss6_3] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

% perturbation to flux 4
exp_pid = 10; % 'V4max'
exp_pval = [0;.1;.3;.5;.7;.9;1];
[xss6_4,fss6_4] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

%% rerun all perturbations with new parameters from x_opt_4
% check if new parameters give the same perturbed flux value
opts.x0 = xss1;
opts.odep = odep_bkp;
opts.tspan = 0:.1:10000;
allopt_id_4 = cat(2,x_opt_4(:).opt_id);
allx_opt_4 = cat(1,x_opt_4(:).x_opt);
opts.odep(allopt_id_4) = allx_opt_4;

% perturbation to flux 1
exp_pid = 11; % 'k1cat'
exp_pval = [.1;.5;1.0;1.5;2];
[xss7_1,fss7_1] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

% perturbation to flux 2
exp_pid = 13; % 'V2max'
exp_pval = [.1;.5;1.0;1.5;2];
[xss7_2,fss7_2] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

% perturbation to flux 3
% restore from backup
exp_pid = 12; % 'V3max'
exp_pval = [.1;.5;1.0;1.5;2];
[xss7_3,fss7_3] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

% perturbation to flux 4
exp_pid = 10; % 'V4max'
exp_pval = [0;.1;.3;.5;.7;.9;1];
[xss7_4,fss7_4] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

% [~,~,xss4_1,fss4_1] = check_kin_kotte(opts);

%% compare all optimized parameter perturbations with original model 
% perturbations
% xss_1,xss_2,xss_3,xss_4 - perturbed fluxes
fprintf('\n==================================================\n');
fprintf('Result Diagnostics\n');
fprintf('==================================================\n');
fprintf('Data Set\tPerturbation\tL2-norm\n');
fprintf('==================================================\n');
dataset = {'x_opt_1','x_opt_2','x_opt_3','x_opt_4'};
perturbations = {'flux 1','flux 2','flux 3','flux 4'};

% x_opt_1
% parameter norms
p_1_diff = norm(allx_opt_1-odep_bkp(allopt_id_1)');

fss4_1_diff = norm(fss4_1-fss_1);
fss4_2_diff = norm(fss4_2-fss_2);
fss4_3_diff = norm(fss4_3-fss_3);
fss4_4_diff = norm(fss4_4-fss_4);
fprintf('%s \t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n',...
         dataset{1},perturbations{1},fss4_1_diff,...
                    perturbations{2},fss4_2_diff,...
                    perturbations{3},fss4_3_diff,...
                    perturbations{4},fss4_4_diff);
fprintf('Parameter Norm : %6.4g\n',p_1_diff);                

% x_opt_2
% parameter norms
p_2_diff = norm(allx_opt_2-odep_bkp(allopt_id_2)');

fss5_1_diff = norm(fss5_1-fss_1);
fss5_2_diff = norm(fss5_2-fss_2);
fss5_3_diff = norm(fss5_3-fss_3);
fss5_4_diff = norm(fss5_4-fss_4);
fprintf('%s \t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n',...
         dataset{2},perturbations{1},fss5_1_diff,...
                    perturbations{2},fss5_2_diff,...
                    perturbations{3},fss5_3_diff,...
                    perturbations{4},fss5_4_diff);
fprintf('Parameter Norm : %6.4g\n',p_2_diff);

% x_opt_3
% parameter norms
p_3_diff = norm(allx_opt_3-odep_bkp(allopt_id_3)');

fss6_1_diff = norm(fss6_1-fss_1);
fss6_2_diff = norm(fss6_2-fss_2);
fss6_3_diff = norm(fss6_3-fss_3);
fss6_4_diff = norm(fss6_4-fss_4);
fprintf('%s \t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n',...
         dataset{3},perturbations{1},fss6_1_diff,...
                    perturbations{2},fss6_2_diff,...
                    perturbations{3},fss6_3_diff,...
                    perturbations{4},fss6_4_diff);
fprintf('Parameter Norm : %6.4g\n',p_3_diff);

% x_opt_4
% parameter norms
p_4_diff = norm(allx_opt_4-odep_bkp(allopt_id_4)');

fss7_1_diff = norm(fss7_1-fss_1);
fss7_2_diff = norm(fss7_2-fss_2);
fss7_3_diff = norm(fss7_3-fss_3);
fss7_4_diff = norm(fss7_4-fss_4);
fprintf('%s \t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n',...
         dataset{4},perturbations{1},fss7_1_diff,...
                    perturbations{2},fss7_2_diff,...
                    perturbations{3},fss7_3_diff,...
                    perturbations{4},fss7_4_diff);  
fprintf('Parameter Norm : %6.4g\n',p_4_diff);
fprintf('================================================\n');                
