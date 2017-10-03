% identifiability analysis with ss data for k1cat
%% load mle data
if ~exist('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel')
    status = 2;
    fprintf('\nLinux System\n');
else 
    status = 1;
    fprintf('\nWindows System\n');
end

if status == 1    
    load('C:/Users/shyam/Documents/Courses/CHE1125Project/results/ple/mle_all_acetate_oct2.mat');
elseif status == 2    
    load('~/Documents/Courses/CHE1125Project/results/ple/mle_all_acetate_oct2.mat');    
end

% use wt as initial value for all perturbations
xinit = noisy_xss(:,1);
yinit = noisy_fss(:,1);
input_data = exp_select_sol.exp_pval;

freq = .1;
ynoise_var = .01;

optim_opts = struct('pname','k1cat','nc',3,'nf',6,'npert',1,...                    
                    'plim',[0.001 6],...
                    'minmax_step',[1e-6 .5],...
                    'casmodelfun',@kotteCASident_pert,...
                    'xinit',xinit,...
                    'yinit',yinit,...
                    'xexp',exp_select_sol.xss,...
                    'yexp',exp_select_sol.fss([1 3 4 5],:),...
                    'ynoise_var',ynoise_var,...
                    'input_data',input_data,...
                    'integratorfun','RK4integrator_cas',...
                    'odep',odep_bkp,...
                    'tspan',opts.tspan,...
                    'freq',freq,...
                    'x0',noisy_xss(:,1));

% set confidence interval threshold for PLE 
alpha = .90; % alpha quantile for chi2 distribution
dof = 12; % degrees of freedom
% chi2 alpha quantile
delta_alpha_1 = chi2inv(alpha,1);      
delta_alpha_all = chi2inv(alpha,dof);
thetai_fixed_value = MLE_noisy.mle_pval(11);
theta_step = 0;

% loop all the abopve statements for complete identifiability algforithm
maxiter = 1000;

% initial value for optimization
scale = ones(8,1);

% p0 for k1cat
scale(3) = 1e6;
p0 = [opts.odep(1:5)';opts.odep(10);opts.odep(12:13)']./scale;

%% call PLE evaluation function
pos_neg = [1 3];
nid = length(pos_neg);
PLEvals = cell(nid,1);
for id = 1:1 % nid
    PLEvals{id} =...
    getPLE(thetai_fixed_value,theta_step,p0,opts.odep,...
           delta_alpha_1,optim_opts,maxiter,pos_neg(id));             
end

%% collect data and plot from parallel estimation
PLE_unify = unifyPLEres(PLEvals);
plotPLE(PLE_unify,delta_alpha_1,delta_alpha_all);  
