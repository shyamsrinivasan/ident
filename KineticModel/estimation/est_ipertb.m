% estimate parameters for same model (kotte) after perturbations
% parameters are estimated on the basis of a single set of perturbations
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

%% use only perturbation set to flux 1 to get parameters
x_opt_1 = optimize_p(opts,xss_1,fss_1,plist,odep_bkp);

%% use only perturbation set to flux 2 to get parameters
x_opt_5 = optimize_p(opts,xss_2,fss_2,plist,odep_bkp);

%% use only perturbation set to flux 3 to get parameters
x_opt_6 = optimize_p(opts,xss_3,fss_3,plist,odep_bkp);

%% use only perturbation set to flux 4 to get parameters
x_opt_7 = optimize_p(opts,xss_4,fss_4,plist,odep_bkp);

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

%% rerun all perturbations with new parameters from x_opt_5
% check if new parameters give the same perturbed flux value
opts.x0 = xss1;
opts.odep = odep_bkp;
opts.tspan = 0:.1:10000;
allopt_id_5 = cat(2,x_opt_5(:).opt_id);
allx_opt_5 = cat(1,x_opt_5(:).x_opt);
opts.odep(allopt_id_5) = allx_opt_5;

% perturbation to flux 1
exp_pid = 11; % 'k1cat'
exp_pval = [.1;.5;1.0;1.5;2];
[xss9_1,fss9_1] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

% perturbation to flux 2
exp_pid = 13; % 'V2max'
exp_pval = [.1;.5;1.0;1.5;2];
[xss9_2,fss9_2] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

% perturbation to flux 3
% restore from backup
exp_pid = 12; % 'V3max'
exp_pval = [.1;.5;1.0;1.5;2];
[xss9_3,fss9_3] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

% perturbation to flux 4
exp_pid = 10; % 'V4max'
exp_pval = [0;.1;.3;.5;.7;.9;1];
[xss9_4,fss9_4] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

%% rerun all perturbations with new parameters from x_opt_6
% check if new parameters give the same perturbed flux value
opts.x0 = xss1;
opts.odep = odep_bkp;
opts.tspan = 0:.1:10000;
allopt_id_6 = cat(2,x_opt_6(:).opt_id);
allx_opt_6 = cat(1,x_opt_6(:).x_opt);
opts.odep(allopt_id_6) = allx_opt_6;

% perturbation to flux 1
exp_pid = 11; % 'k1cat'
exp_pval = [.1;.5;1.0;1.5;2];
[xss10_1,fss10_1] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

% perturbation to flux 2
exp_pid = 13; % 'V2max'
exp_pval = [.1;.5;1.0;1.5;2];
[xss10_2,fss10_2] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

% perturbation to flux 3
% restore from backup
exp_pid = 12; % 'V3max'
exp_pval = [.1;.5;1.0;1.5;2];
[xss10_3,fss10_3] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

% perturbation to flux 4
exp_pid = 10; % 'V4max'
exp_pval = [0;.1;.3;.5;.7;.9;1];
[xss10_4,fss10_4] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

%% rerun all perturbations with new parameters from x_opt_7
% check if new parameters give the same perturbed flux value
opts.x0 = xss1;
opts.odep = odep_bkp;
opts.tspan = 0:.1:10000;
allopt_id_7 = cat(2,x_opt_7(:).opt_id);
allx_opt_7 = cat(1,x_opt_7(:).x_opt);
opts.odep(allopt_id_7) = allx_opt_7;

% perturbation to flux 1
exp_pid = 11; % 'k1cat'
exp_pval = [.1;.5;1.0;1.5;2];
[xss11_1,fss11_1] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

% perturbation to flux 2
exp_pid = 13; % 'V2max'
exp_pval = [.1;.5;1.0;1.5;2];
[xss11_2,fss11_2] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

% perturbation to flux 3
% restore from backup
exp_pid = 12; % 'V3max'
exp_pval = [.1;.5;1.0;1.5;2];
[xss11_3,fss11_3] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

% perturbation to flux 4
exp_pid = 10; % 'V4max'
exp_pval = [0;.1;.3;.5;.7;.9;1];
[xss11_4,fss11_4] = runperturbations(@perturb_nonoise,exp_pid,exp_pval,opts);
close all

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
p_5_diff = norm(allx_opt_5-odep_bkp(allopt_id_5)');

fss9_1_diff = norm(fss9_1-fss_1);
fss9_2_diff = norm(fss9_2-fss_2);
fss9_3_diff = norm(fss9_3-fss_3);
fss9_4_diff = norm(fss9_4-fss_4);
fprintf('%s \t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n',...
         dataset{2},perturbations{1},fss9_1_diff,...
                    perturbations{2},fss9_2_diff,...
                    perturbations{3},fss9_3_diff,...
                    perturbations{4},fss9_4_diff);
fprintf('Parameter Norm : %6.4g\n',p_5_diff);

% x_opt_3
% parameter norms
p_6_diff = norm(allx_opt_6-odep_bkp(allopt_id_6)');

fss10_1_diff = norm(fss10_1-fss_1);
fss10_2_diff = norm(fss10_2-fss_2);
fss10_3_diff = norm(fss10_3-fss_3);
fss10_4_diff = norm(fss10_4-fss_4);
fprintf('%s \t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n',...
         dataset{3},perturbations{1},fss10_1_diff,...
                    perturbations{2},fss10_2_diff,...
                    perturbations{3},fss10_3_diff,...
                    perturbations{4},fss10_4_diff);
fprintf('Parameter Norm : %6.4g\n',p_6_diff);

% x_opt_4
% parameter norms
p_7_diff = norm(allx_opt_7-odep_bkp(allopt_id_7)');

fss11_1_diff = norm(fss11_1-fss_1);
fss11_2_diff = norm(fss11_2-fss_2);
fss11_3_diff = norm(fss11_3-fss_3);
fss11_4_diff = norm(fss11_4-fss_4);
fprintf('%s \t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n',...
         dataset{4},perturbations{1},fss11_1_diff,...
                    perturbations{2},fss11_2_diff,...
                    perturbations{3},fss11_3_diff,...
                    perturbations{4},fss11_4_diff);  
fprintf('Parameter Norm : %6.4g\n',p_7_diff);
fprintf('================================================\n');            