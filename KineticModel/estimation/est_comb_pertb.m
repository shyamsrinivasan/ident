% try different combinations of perturbations to predict enzyme kinetic
% parameters
plist = {'K1ac','K3fdp','L3fdp','K3pep','K2pep','vemax','KeFDP','ne',...
        'd','V4max','k1cat','V3max','V2max','K1pep','K2fdp','rhoA'};

% generate experimental data - get initial ss
tspan = 0:0.1:300;
[xdyn,fdyn,xss1,fss1,opts] = run_nonoise(tspan);

% backup parameters and initial conditions
ival_bkp = opts.x0;
odep_bkp = opts.odep;

%% perturbation to all fluxes 
opts.x0 = xss1;
opts.tspan = 0:.1:300;
opts.odep = odep_bkp;
% flux 1, 2 and 3 % k1cat, 'V2max', 'V3max'
ptopts = struct('exp_pid',{11,13,12},...
                'exp_pval',{[.5;1.0;1.5;2],...
                            [.1;.5;1.0;1.5;2],...
                            [.1;.5;1.0;1.5;2]}); 
sol = getperturbations(ptopts,@perturb_nonoise,opts);
close all

% flux 4 V4max
opts.tspan = 0:.1:10000;
ptopts = struct('exp_pid',10,'exp_pval',[0;.1;.3;.5;.7;.9;1]);
sol = getperturbations(ptopts,@perturb_nonoise,opts,sol);
close all

%% use different combinations of perturbation sets to get parameters
optimopts = struct('xss',{cat(2,sol(1).xss,sol(3).xss),...
                          cat(2,sol(2).xss,sol(4).xss),...
                          cat(2,sol(1).xss,sol(4).xss),...
                          cat(2,sol(2).xss,sol(3).xss),...
                          cat(2,sol(3).xss,sol(4).xss),...
                          cat(2,sol(1).xss,sol(2).xss,sol(4).xss),...
                          cat(2,sol(2).xss,sol(3).xss,sol(4).xss)},...
                   'fss',{cat(2,sol(1).fss,sol(3).fss),...
                          cat(2,sol(2).fss,sol(4).fss),...
                          cat(2,sol(1).fss,sol(4).fss),...
                          cat(2,sol(2).fss,sol(3).fss),...
                          cat(2,sol(3).fss,sol(4).fss),...
                          cat(2,sol(1).fss,sol(2).fss,sol(4).fss),...
                          cat(2,sol(2).fss,sol(3).fss,sol(4).fss)});
opt_sol = runoptimp(opts,plist,odep_bkp,optimopts);    

%% rerun all perturbations with new parameters from opt_sol
% check if new parameters give the same perturbed flux value
opts.x0 = xss1;
opts.odep = odep_bkp;
opts.tspan = 0:.1:10000;
allopt_id_8 = cat(2,opt_sol{1}(:).opt_id);
allx_opt_8 = cat(1,opt_sol{1}(:).x_opt);
opts.odep(allopt_id_8) = allx_opt_8;

% flux 1, 2 and 3 % k1cat, 'V2max', 'V3max'
ptopts = struct('exp_pid',{11,13,12},...
                'exp_pval',{[.5;1.0;1.5;2],...
                            [.1;.5;1.0;1.5;2],...
                            [.1;.5;1.0;1.5;2]}); 
sol_opt_p_8 = getperturbations(ptopts,@perturb_nonoise,opts);
ptopts = struct('exp_pid',10,'exp_pval',[0;.1;.3;.5;.7;.9;1]);
sol_opt_p_8 = getperturbations(ptopts,@perturb_nonoise,opts,sol_opt_p_8);
close all

%% rerun all perturbations with new parameters from x_opt_2
% check if new parameters give the same perturbed flux value
opts.x0 = xss1;
opts.odep = odep_bkp;
opts.tspan = 0:.1:10000;
allopt_id_9 = cat(2,opt_sol{2}(:).opt_id);
allx_opt_9 = cat(1,opt_sol{2}(:).x_opt);
opts.odep(allopt_id_9) = allx_opt_9;

% flux 1, 2 and 3 % k1cat, 'V2max', 'V3max'
ptopts = struct('exp_pid',{11,13,12},...
                'exp_pval',{[.5;1.0;1.5;2],...
                            [.1;.5;1.0;1.5;2],...
                            [.1;.5;1.0;1.5;2]}); 
sol_opt_p_9 = getperturbations(ptopts,@perturb_nonoise,opts);
ptopts = struct('exp_pid',10,'exp_pval',[0;.1;.3;.5;.7;.9;1]);
sol_opt_p_9 = getperturbations(ptopts,@perturb_nonoise,opts,sol_opt_p_9);
close all

%% rerun all perturbations with new parameters from x_opt_3
% check if new parameters give the same perturbed flux value
opts.x0 = xss1;
opts.odep = odep_bkp;
opts.tspan = 0:.1:10000;
allopt_id_10 = cat(2,opt_sol{3}(:).opt_id);
allx_opt_10 = cat(1,opt_sol{3}(:).x_opt);
opts.odep(allopt_id_10) = allx_opt_10;

% flux 1, 2 and 3 % k1cat, 'V2max', 'V3max'
ptopts = struct('exp_pid',{11,13,12},...
                'exp_pval',{[.1;.5;1.0;1.5;2],...
                            [.1;.5;1.0;1.5;2],...
                            [.1;.5;1.0;1.5;2]}); 
sol_opt_p_10 = getperturbations(ptopts,@perturb_nonoise,opts);
ptopts = struct('exp_pid',10,'exp_pval',[0;.1;.3;.5;.7;.9;1]);
sol_opt_p_10 = getperturbations(ptopts,@perturb_nonoise,opts,sol_opt_p_10);
close all

%% rerun all perturbations with new parameters from x_opt_4
% check if new parameters give the same perturbed flux value
opts.x0 = xss1;
opts.odep = odep_bkp;
opts.tspan = 0:.1:10000;
allopt_id_11 = cat(2,opt_sol{4}(:).opt_id);
allx_opt_11 = cat(1,opt_sol{4}(:).x_opt);
opts.odep(allopt_id_11) = allx_opt_11;

% flux 1, 2 and 3 % k1cat, 'V2max', 'V3max'
ptopts = struct('exp_pid',{11,13,12},...
                'exp_pval',{[.1;.5;1.0;1.5;2],...
                            [.1;.5;1.0;1.5;2],...
                            [.1;.5;1.0;1.5;2]}); 
sol_opt_p_11 = getperturbations(ptopts,@perturb_nonoise,opts);
ptopts = struct('exp_pid',10,'exp_pval',[0;.1;.3;.5;.7;.9;1]);
sol_opt_p_11 = getperturbations(ptopts,@perturb_nonoise,opts,sol_opt_p_11);
close all

%% rerun all perturbations with new parameters from x_opt_4
% check if new parameters give the same perturbed flux value
opts.x0 = xss1;
opts.odep = odep_bkp;
opts.tspan = 0:.1:10000;
allopt_id_12 = cat(2,opt_sol{5}(:).opt_id);
allx_opt_12 = cat(1,opt_sol{5}(:).x_opt);
opts.odep(allopt_id_12) = allx_opt_12;

% flux 1, 2 and 3 % k1cat, 'V2max', 'V3max'
ptopts = struct('exp_pid',{11,13,12},...
                'exp_pval',{[.1;.5;1.0;1.5;2],...
                            [.1;.5;1.0;1.5;2],...
                            [.1;.5;1.0;1.5;2]}); 
sol_opt_p_12 = getperturbations(ptopts,@perturb_nonoise,opts);
ptopts = struct('exp_pid',10,'exp_pval',[0;.1;.3;.5;.7;.9;1]);
sol_opt_p_12 = getperturbations(ptopts,@perturb_nonoise,opts,sol_opt_p_12);
close all

%% rerun all perturbations with new parameters from x_opt_4
% check if new parameters give the same perturbed flux value
opts.x0 = xss1;
opts.odep = odep_bkp;
opts.tspan = 0:.1:10000;
allopt_id_13 = cat(2,opt_sol{6}(:).opt_id);
allx_opt_13 = cat(1,opt_sol{6}(:).x_opt);
opts.odep(allopt_id_13) = allx_opt_13;

% flux 1, 2 and 3 % k1cat, 'V2max', 'V3max'
ptopts = struct('exp_pid',{11,13,12},...
                'exp_pval',{[.1;.5;1.0;1.5;2],...
                            [.1;.5;1.0;1.5;2],...
                            [.1;.5;1.0;1.5;2]}); 
sol_opt_p_13 = getperturbations(ptopts,@perturb_nonoise,opts);
ptopts = struct('exp_pid',10,'exp_pval',[0;.1;.3;.5;.7;.9;1]);
sol_opt_p_13 = getperturbations(ptopts,@perturb_nonoise,opts,sol_opt_p_13);
close all

%% rerun all perturbations with new parameters from x_opt_4
% check if new parameters give the same perturbed flux value
opts.x0 = xss1;
opts.odep = odep_bkp;
opts.tspan = 0:.1:10000;
allopt_id_14 = cat(2,opt_sol{7}(:).opt_id);
allx_opt_14 = cat(1,opt_sol{7}(:).x_opt);
opts.odep(allopt_id_14) = allx_opt_14;

% flux 1, 2 and 3 % k1cat, 'V2max', 'V3max'
ptopts = struct('exp_pid',{11,13,12},...
                'exp_pval',{[.1;.5;1.0;1.5;2],...
                            [.1;.5;1.0;1.5;2],...
                            [.1;.5;1.0;1.5;2]}); 
sol_opt_p_14 = getperturbations(ptopts,@perturb_nonoise,opts);
ptopts = struct('exp_pid',10,'exp_pval',[0;.1;.3;.5;.7;.9;1]);
sol_opt_p_14 = getperturbations(ptopts,@perturb_nonoise,opts,sol_opt_p_14);
close all

%% compare all optimized parameter perturbations with original model 
% perturbations
% xss_1,xss_2,xss_3,xss_4 - perturbed fluxes
fprintf('\n==================================================\n');
fprintf('Result Diagnostics\n');
fprintf('==================================================\n');
fprintf('Data Set\tPerturbation\tL2-norm\n');
fprintf('==================================================\n');
dataset = {'x_opt_1','x_opt_2','x_opt_3','x_opt_4','x_opt_5','x_opt_6','x_opt_7'};
perturbations = {'flux 1','flux 2','flux 3','flux 4'};

% x_opt_1
% parameter norms
p_8_diff = norm(allx_opt_8-odep_bkp(allopt_id_8)');

% fss4_1_diff = norm(sol_opt_p_1(1).fss-sol(1).fss);
% fss4_2_diff = norm(sol_opt_p_1(2).fss-sol(2).fss);
% fss4_3_diff = norm(sol_opt_p_1(3).fss-sol(3).fss);
% fss4_4_diff = 1.0; % norm(fss4_4-fss_4);
fprintf('%s \t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n',...
         dataset{1},perturbations{1},norm(sol_opt_p_8(1).fss-sol(1).fss),...
                    perturbations{2},norm(sol_opt_p_8(2).fss-sol(2).fss),...
                    perturbations{3},norm(sol_opt_p_8(3).fss-sol(3).fss),...
                    perturbations{4},norm(sol_opt_p_8(4).fss-sol(4).fss));
fprintf('Parameter Norm : %6.4g\n',p_8_diff);                

% x_opt_2
% parameter norms
p_9_diff = norm(allx_opt_9-odep_bkp(allopt_id_9)');

% fss5_1_diff = norm(sol_opt_p_2(1).fss-sol(1).fss);
% fss5_2_diff = norm(sol_opt_p_2(2).fss-sol(2).fss);
% fss5_3_diff = norm(sol_opt_p_2(3).fss-sol(3).fss);
% fss5_4_diff = 1.0; % norm(fss5_4-fss_4);
fprintf('%s \t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n',...
         dataset{2},perturbations{1},norm(sol_opt_p_9(1).fss-sol(1).fss),...
                    perturbations{2},norm(sol_opt_p_9(2).fss-sol(2).fss),...
                    perturbations{3},norm(sol_opt_p_9(3).fss-sol(3).fss),...
                    perturbations{4},norm(sol_opt_p_9(4).fss-sol(4).fss));
fprintf('Parameter Norm : %6.4g\n',p_9_diff);

% x_opt_3
% parameter norms
p_10_diff = norm(allx_opt_10-odep_bkp(allopt_id_10)');

% fss6_1_diff = norm(sol_opt_p_3(1).fss-sol(1).fss);
% fss6_2_diff = norm(sol_opt_p_3(2).fss-sol(2).fss);
% fss6_3_diff = norm(sol_opt_p_3(3).fss-sol(3).fss);
% fss6_4_diff = 1.0; % norm(fss6_4-fss_4);
fprintf('%s \t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n',...
         dataset{3},perturbations{1},norm(sol_opt_p_10(1).fss-sol(1).fss),...
                    perturbations{2},norm(sol_opt_p_10(2).fss-sol(2).fss),...
                    perturbations{3},norm(sol_opt_p_10(3).fss-sol(3).fss),...
                    perturbations{4},norm(sol_opt_p_10(4).fss-sol(4).fss));
fprintf('Parameter Norm : %6.4g\n',p_10_diff);

% x_opt_4
% parameter norms
p_11_diff = norm(allx_opt_11-odep_bkp(allopt_id_11)');

% fss7_1_diff = norm(sol_opt_p_4(1).fss-sol(1).fss);
% fss7_2_diff = norm(sol_opt_p_4(2).fss-sol(2).fss);
% fss7_3_diff = norm(sol_opt_p_4(3).fss-sol(3).fss);
% fss7_4_diff = 1.0; % norm(fss7_4-fss_4);
fprintf('%s \t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n',...
         dataset{4},perturbations{1},norm(sol_opt_p_11(1).fss-sol(1).fss),...
                    perturbations{2},norm(sol_opt_p_11(2).fss-sol(2).fss),...
                    perturbations{3},norm(sol_opt_p_11(3).fss-sol(3).fss),...
                    perturbations{4},norm(sol_opt_p_11(4).fss-sol(4).fss));  
fprintf('Parameter Norm : %6.4g\n',p_11_diff);
fprintf('================================================\n');   

% x_opt_4
% parameter norms
p_12_diff = norm(allx_opt_12-odep_bkp(allopt_id_12)');

% fss7_1_diff = norm(sol_opt_p_4(1).fss-sol(1).fss);
% fss7_2_diff = norm(sol_opt_p_4(2).fss-sol(2).fss);
% fss7_3_diff = norm(sol_opt_p_4(3).fss-sol(3).fss);
% fss7_4_diff = 1.0; % norm(fss7_4-fss_4);
fprintf('%s \t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n',...
         dataset{5},perturbations{1},norm(sol_opt_p_12(1).fss-sol(1).fss),...
                    perturbations{2},norm(sol_opt_p_12(2).fss-sol(2).fss),...
                    perturbations{3},norm(sol_opt_p_12(3).fss-sol(3).fss),...
                    perturbations{4},norm(sol_opt_p_12(4).fss-sol(4).fss));  
fprintf('Parameter Norm : %6.4g\n',p_12_diff);
fprintf('================================================\n');   

% x_opt_4
% parameter norms
p_13_diff = norm(allx_opt_13-odep_bkp(allopt_id_13)');

% fss7_1_diff = norm(sol_opt_p_4(1).fss-sol(1).fss);
% fss7_2_diff = norm(sol_opt_p_4(2).fss-sol(2).fss);
% fss7_3_diff = norm(sol_opt_p_4(3).fss-sol(3).fss);
% fss7_4_diff = 1.0; % norm(fss7_4-fss_4);
fprintf('%s \t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n',...
         dataset{6},perturbations{1},norm(sol_opt_p_13(1).fss-sol(1).fss),...
                    perturbations{2},norm(sol_opt_p_13(2).fss-sol(2).fss),...
                    perturbations{3},norm(sol_opt_p_13(3).fss-sol(3).fss),...
                    perturbations{4},norm(sol_opt_p_13(4).fss-sol(4).fss));  
fprintf('Parameter Norm : %6.4g\n',p_13_diff);
fprintf('================================================\n');   

% x_opt_4
% parameter norms
p_14_diff = norm(allx_opt_14-odep_bkp(allopt_id_14)');

% fss7_1_diff = norm(sol_opt_p_4(1).fss-sol(1).fss);
% fss7_2_diff = norm(sol_opt_p_4(2).fss-sol(2).fss);
% fss7_3_diff = norm(sol_opt_p_4(3).fss-sol(3).fss);
% fss7_4_diff = 1.0; % norm(fss7_4-fss_4);
fprintf('%s \t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n',...
         dataset{7},perturbations{1},norm(sol_opt_p_14(1).fss-sol(1).fss),...
                    perturbations{2},norm(sol_opt_p_14(2).fss-sol(2).fss),...
                    perturbations{3},norm(sol_opt_p_14(3).fss-sol(3).fss),...
                    perturbations{4},norm(sol_opt_p_14(4).fss-sol(4).fss));  
fprintf('Parameter Norm : %6.4g\n',p_14_diff);
fprintf('================================================\n');   


