% generate experimental data - get initial ss
tspan = 0:.1:300;
[xdyn,fdyn,xss,fss,opts] = run_nonoise(tspan);

% backup parameters and initial conditions
ival_bkp = opts.x0;
odep_bkp = opts.odep;

% perturb system from non-noisy initial conditions
opts.tspan = 0:.1:300;
pt_sol_id = [1 2 3];
[exp_sol,no_noise_sol] = dopert_nonoise_dyn(opts,xss,fss,odep_bkp,pt_sol_id);
close all
% load('C:/Users/shyam/Documents/Courses/CHE1125Project/IntegratedModels/KineticModel/estimation/noiseless_model/pdata_nn_aug29');

freq = 1:200:3001;
optim_opts = struct('pname','k1cat','nc',3,'nf',6,'npert',2,...
                    'plim',[0.001 6],...
                    'casmodelfun',@kotteCASident_pert,...
                    'integratorfun','RK4integrator_cas',...
                    'odep',odep_bkp,...
                    'tspan',opts.tspan,...
                    'freq',freq,...
                    'x0',xss,...
                    'xinit',exp_sol.xdyn(:,1),...
                    'xexp',exp_sol.xdyn(:,freq));

% set confidence interval threshold for PLE 
alpha = .90; % alpha quantile for chi2 distribution
dof = 1; % degrees of freedom
% chi2 alpha quantile
delta_alpha = chi2inv(alpha,dof);      
thetai_fixed_value = .1;
theta_step = 0;

% loop all the abopve statements for complete identifiability algforithm
maxiter = 20;

% initial value for optimization
scale = ones(12,1);

% p0 for K1ac
% scale(2) = 1e6;
% p0 = opts.odep(2:13)'./scale;

% p0 for k1cat
scale(3) = 1e6;
p0 = [opts.odep(1:10)';opts.odep(12:13)']./scale;

% call PLE evaluation function
[PLEvals] =...
getPLE(thetai_fixed_value,theta_step,p0,opts.odep,delta_alpha,optim_opts,maxiter,2);

% get confidence interval bounds from pl data
% sigma = pleCI_finite_sample(delta_alpha,PLEvals_pos.chiPLE,PLEvals_pos.thetai_inc);

% figure
% hold on

% figure for PLE
%     line(thetai_inc,chiPLE(iter),'LineStyle','none','Marker','.','MarkerSize',10);

