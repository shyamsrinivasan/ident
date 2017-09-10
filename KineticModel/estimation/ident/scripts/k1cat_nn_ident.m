% identifiability analysis with ss data for k1cat

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
use_pert = 1;
npert = length(use_pert);
[exp_select_sol,no_noise_select_sol] = parseperturbations(no_noise_sol,use_pert);

% use wt as initial value for all perturbations
xinit = repmat(xss,npert,1);
yinit = repmat(fss,npert,1);

freq = [1:50:1500 1501:1500:3001];
optim_opts = struct('pname','k1cat','nc',3,'nf',6,'npert',npert,...                    
                    'plim',[0.001 6],...
                    'minmax_step',[1e-6 .4],...
                    'casmodelfun',@kotteCASident_pert,...
                    'integratorfun','RK4integrator_cas',...
                    'odep',odep_bkp,...
                    'tspan',opts.tspan,...
                    'freq',freq,...
                    'x0',xss,...
                    'xinit',xinit,...
                    'yinit',yinit,...
                    'xexp',exp_select_sol.xdyn(:,freq),...
                    'yexp',exp_select_sol.fdyn([1 3 4 5],freq));

% set confidence interval threshold for PLE 
alpha = .90; % alpha quantile for chi2 distribution
dof = 12; % degrees of freedom
% chi2 alpha quantile
delta_alpha_1 = chi2inv(alpha,1);      
delta_alpha_all = chi2inv(alpha,dof);
thetai_fixed_value = .1;
theta_step = 0;

% loop all the abopve statements for complete identifiability algforithm
maxiter = 1000;

% initial value for optimization
scale = ones(12,1);

% p0 for K1ac
% scale(2) = 1e6;
% p0 = opts.odep(2:13)'./scale;

% p0 for k1cat
scale(3) = 1e6;
p0 = [opts.odep(1:10)';opts.odep(12:13)']./scale;

%% call PLE evaluation function
pos_neg = [1 3];
nid = length(pos_neg);
PLEvals = cell(nid,1);
parfor id = 1:nid
    PLEvals{id} =...
    getPLE(thetai_fixed_value,theta_step,p0,opts.odep,...
           delta_alpha_1,optim_opts,maxiter,pos_neg(id));

    plotPLE(PLEvals{id},delta_alpha_1,delta_alpha_all);                
end
%%    



