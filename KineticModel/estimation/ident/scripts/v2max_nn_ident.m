% identifiability analysis with ss data for V2max

%% load noise free data
if ~exist('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel')
    status = 2;
    fprintf('\nLinux System\n');
else 
    status = 1;
    fprintf('\nWindows System\n');
end

if status == 1    
    load('C:/Users/shyam/Documents/Courses/CHE1125Project/IntegratedModels/KineticModel/estimation/noiseless_model/pdata_nn_sep1');
elseif status == 2    
    load('~/Documents/Courses/CHE1125Project/IntegratedModels/KineticModel/estimation/noiseless_model/pdata_nn_sep1');    
end

%% set PLE options
% collect only needed perturbations for analysis
avail_pert = size(no_noise_sol,2);
use_pert = [1 2 3 avail_pert];
npert = length(use_pert);
[exp_select_sol,no_noise_select_sol] = parseperturbations(no_noise_sol,use_pert);

% use wt as initial value for all perturbations
xinit = repmat(exp_select_sol.xss(:,end),npert,1);

freq = 1:3000:3001;
optim_opts = struct('pname','V2max','nc',3,'nf',6,'npert',npert,...
                    'nunpert',10,...
                    'plim',[0.001 6],...
                    'casmodelfun',@kotteCASident_pert,...
                    'integratorfun','RK4integrator_cas',...
                    'odep',odep_bkp,...
                    'tspan',opts.tspan,...
                    'freq',freq,...
                    'x0',xss,...
                    'xinit',xinit,...
                    'xexp',exp_select_sol.xdyn(:,freq),...
                    'p_pert',exp_select_sol.p_pert);

% set confidence interval threshold for PLE 
alpha = .90; % alpha quantile for chi2 distribution
dof = 1; % degrees of freedom
% chi2 alpha quantile
delta_alpha = chi2inv(alpha,dof);      
thetai_fixed_value = .1;
theta_step = 0;

% loop all the abopve statements for complete identifiability algforithm
maxiter = 1000;

% initial value for optimization
scale = ones(10,1);

% p0 for K1ac
% scale(2) = 1e6;
% p0 = opts.odep(2:13)'./scale;

% p0 for k1cat
scale(3) = 1e6;
p0 = [opts.odep(1:10)'./scale;repmat(opts.odep(11:12)',npert,1)];
identobj(p0,.1,optim_opts);

%% call PLE evaluation function
[PLEvals] =...
getPLE(thetai_fixed_value,theta_step,p0,opts.odep,delta_alpha,optim_opts,maxiter,2);
                
%%    



