% estimate parameters for same model (kotte) after perturbations
% parameters are estimated on the basis of a single set of perturbations
% all parameters available
plist = {'K1ac','K3fdp','L3fdp','K3pep','K2pep','vemax','KeFDP','ne',...
        'd','V4max','k1cat','V3max','V2max','K1pep','K2fdp','rhoA','acetate'};

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

%% use single perturbation sets to get parameters
optimopts = struct('xss',{sol(1).xss,sol(2).xss,...
                          sol(3).xss,sol(4).xss},...
                   'fss',{sol(1).fss,sol(2).fss,...
                          sol(3).fss,sol(4).fss});
opt_sol = runoptimp(opts,plist,odep_bkp,optimopts,@optimize_p);                       
% x_opt_1 = optimize_p(opts,xss_1,fss_1,plist,odep_bkp);

%% rerun all perturbations with new parameters from x_opt_1
% check if new parameters give the same perturbed flux value
opts.x0 = xss1;
opts.odep = odep_bkp;
opts.tspan = 0:.1:10000;
allopt_id_1 = cat(2,opt_sol{1}(:).opt_id);
allx_opt_1 = cat(1,opt_sol{1}(:).x_opt);
opts.odep(allopt_id_1) = allx_opt_1;

% perturbations to flux 1, 2 and 3 % k1cat, 'V2max', 'V3max'
ptopts = struct('exp_pid',{11,13,12},...
                'exp_pval',{[.5;1.0;1.5;2],...
                            [.1;.5;1.0;1.5;2],...
                            [.1;.5;1.0;1.5;2]}); 
sol_opt_p_1 = getperturbations(ptopts,@perturb_nonoise,opts);
ptopts = struct('exp_pid',10,'exp_pval',[0;.1;.3;.5;.7;.9;1]);
sol_opt_p_1 = getperturbations(ptopts,@perturb_nonoise,opts,sol_opt_p_1);
close all

%% rerun all perturbations with new parameters from x_opt_5
% check if new parameters give the same perturbed flux value
opts.x0 = xss1;
opts.odep = odep_bkp;
opts.tspan = 0:.1:10000;
allopt_id_5 = cat(2,opt_sol{2}(:).opt_id);
allx_opt_5 = cat(1,opt_sol{2}(:).x_opt);
opts.odep(allopt_id_5) = allx_opt_5;

% perturbations to flux 1, 2 and 3 and 4 % k1cat, 'V2max', 'V3max', V4max
ptopts = struct('exp_pid',{11,13,12},...
                'exp_pval',{[.5;1.0;1.5;2],...
                            [.1;.5;1.0;1.5;2],...
                            [.1;.5;1.0;1.5;2]}); 
sol_opt_p_5 = getperturbations(ptopts,@perturb_nonoise,opts);
ptopts = struct('exp_pid',10,'exp_pval',[0;.1;.3;.5;.7;.9;1]);
sol_opt_p_5 = getperturbations(ptopts,@perturb_nonoise,opts,sol_opt_p_5);
close all

%% rerun all perturbations with new parameters from x_opt_6
% check if new parameters give the same perturbed flux value
opts.x0 = xss1;
opts.odep = odep_bkp;
opts.tspan = 0:.1:10000;
allopt_id_6 = cat(2,opt_sol{3}(:).opt_id);
allx_opt_6 = cat(1,opt_sol{3}(:).x_opt);
opts.odep(allopt_id_6) = allx_opt_6;

% flux 1, 2 and 3 and 4 % k1cat, 'V2max', 'V3max'
ptopts = struct('exp_pid',{11,13,12},...
                'exp_pval',{[.5;1.0;1.5;2],...
                            [.1;.5;1.0;1.5;2],...
                            [.1;.5;1.0;1.5;2]}); 
sol_opt_p_6 = getperturbations(ptopts,@perturb_nonoise,opts);
ptopts = struct('exp_pid',10,'exp_pval',[0;.1;.3;.5;.7;.9;1]);
sol_opt_p_6 = getperturbations(ptopts,@perturb_nonoise,opts,sol_opt_p_6);
close all

%% rerun all perturbations with new parameters from x_opt_7
% check if new parameters give the same perturbed flux value
opts.x0 = xss1;
opts.odep = odep_bkp;
opts.tspan = 0:.1:10000;
allopt_id_7 = cat(2,opt_sol{4}(:).opt_id);
allx_opt_7 = cat(1,opt_sol{4}(:).x_opt);
opts.odep(allopt_id_7) = allx_opt_7;

% flux 1, 2 and 3 and 4 % k1cat, 'V2max', 'V3max'
ptopts = struct('exp_pid',{11,13,12},...
                'exp_pval',{[.5;1.0;1.5;2],...
                            [.1;.5;1.0;1.5;2],...
                            [.1;.5;1.0;1.5;2]}); 
sol_opt_p_7 = getperturbations(ptopts,@perturb_nonoise,opts);
ptopts = struct('exp_pid',10,'exp_pval',[0;.1;.3;.5;.7;.9;1]);
sol_opt_p_7 = getperturbations(ptopts,@perturb_nonoise,opts,sol_opt_p_7);
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

fprintf('%s \t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n',...
         dataset{1},perturbations{1},norm(sol_opt_p_1(1).fss-sol(1).fss),...
                    perturbations{2},norm(sol_opt_p_1(2).fss-sol(2).fss),...
                    perturbations{3},norm(sol_opt_p_1(3).fss-sol(3).fss),...
                    perturbations{4},norm(sol_opt_p_1(4).fss-sol(4).fss));
fprintf('Parameter Norm : %6.4g\n',p_1_diff);                    

% x_opt_2
% parameter norms
p_2_diff = norm(allx_opt_5-odep_bkp(allopt_id_5)');

fprintf('%s \t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n',...
         dataset{2},perturbations{1},norm(sol_opt_p_5(1).fss-sol(1).fss),...
                    perturbations{2},norm(sol_opt_p_5(2).fss-sol(2).fss),...
                    perturbations{3},norm(sol_opt_p_5(3).fss-sol(3).fss),...
                    perturbations{4},norm(sol_opt_p_5(4).fss-sol(4).fss));
fprintf('Parameter Norm : %6.4g\n',p_2_diff);

% x_opt_3
% parameter norms
p_3_diff = norm(allx_opt_6-odep_bkp(allopt_id_6)');

fprintf('%s \t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n',...
         dataset{3},perturbations{1},norm(sol_opt_p_6(1).fss-sol(1).fss),...
                    perturbations{2},norm(sol_opt_p_6(2).fss-sol(2).fss),...
                    perturbations{3},norm(sol_opt_p_6(3).fss-sol(3).fss),...
                    perturbations{4},norm(sol_opt_p_6(4).fss-sol(4).fss));
fprintf('Parameter Norm : %6.4g\n',p_3_diff);


% x_opt_4
% parameter norms
p_4_diff = norm(allx_opt_7-odep_bkp(allopt_id_7)');

fprintf('%s \t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n\t\t\t %s \t %6.4g \n',...
         dataset{4},perturbations{1},norm(sol_opt_p_7(1).fss-sol(1).fss),...
                    perturbations{2},norm(sol_opt_p_7(2).fss-sol(2).fss),...
                    perturbations{3},norm(sol_opt_p_7(3).fss-sol(3).fss),...
                    perturbations{4},norm(sol_opt_p_7(4).fss-sol(4).fss));  
fprintf('Parameter Norm : %6.4g\n',p_4_diff);
fprintf('================================================\n');                