% generate experimental data - get initial ss
tspan = 0:.1:300;
[xdyn,fdyn,xss,fss,opts] = run_nonoise(tspan);

% backup parameters and initial conditions
ival_bkp = opts.x0;
odep_bkp = opts.odep;

% perturb system from non-noisy initial conditions
opts.tspan = 0:.1:200;
pt_sol_id = [1 2 3];
[exp_sol,no_noise_sol] = dopert_nonoise_dyn(opts,xss,fss,odep_bkp,pt_sol_id);
close all

optim_opts = struct('pname','K1ac','nc',3,'nf',6,...
                    'casmodelfun',@kotteCASident,...
                    'integratorfun','RK4integrator_cas',...
                    'odep',odep_bkp,...
                    'tspan',opts.tspan,...
                    'x0',xss,...
                    'xinit',no_noise_sol(4).xss,...
                    'xexp',no_noise_sol(1).xdyn(:,1:100:2001));

% set confidence interval threshold for PLE 
alpha = .50; % alpha quantile for chi2 distribution
dof = 1; % degrees of freedom
% chi2 alpha quantile
delta_alpha = chi2inv(alpha,dof);      
thetai_fixed_value = .1;
theta_step = 0;

% loop all the abopve statements for complete identifiability algforithm
maxiter = 100;
chiPLE = zeros(1,maxiter);
xPLE = zeros(12,maxiter);
thetai_inc = zeros(1,maxiter);
thetai_step = zeros(1,maxiter);
obj_step = zeros(1,maxiter);

% initial value for optimization
p0 = opts.odep(2:13)';
iter = 1;

figure
hold on
while iter<maxiter || chiPLE(iter)<delta_alpha
    
    [optsol,thetai_fixed_value,theta_step,obj_new] =...
    PLEiter(thetai_fixed_value,theta_step,p0,opts.odep,delta_alpha,optim_opts);
        
    % new initial p0 = old optimal value
    p0 = optsol.xval;    
     
    % store obj values
    obj_step(iter) = obj_new;
    chiPLE(iter) = optsol.fval;
    xPLE(:,iter) = optsol.xval;
    thetai_inc(iter) = thetai_fixed_value;
    thetai_step(iter) = theta_step;
    
    % figure for PLE
    line(thetai_inc,chiPLE(iter),'LineStyle','none','Marker','.','MarkerSize',10);
    
    iter = iter+1;
end