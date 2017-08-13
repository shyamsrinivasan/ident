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

% problem & objective setup
prob_cas = identopt_setup(optim_opts,.1);

% solve to get optimal parameters
optsol = solve_nlsqopt(prob_cas,opts.odep(2:13)');

% generate new step size for thetai
[theta_step,iter] = adaptive_step(optsol.xval,prob_cas,opts.odep,.1);

% loop all the abopve sttements for complete identifiability algforithm