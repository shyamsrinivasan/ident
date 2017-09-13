% script to generate data (noisy) to be used for analysis 
% in this git repo estimation and identifiability)
% generate experimental data - get initial ss
tspan = 0:0.1:300;
[xdyn,fdyn,xss,fss,opts] = run_nonoise(tspan);

% backup parameters and initial conditions
ival_bkp = opts.x0;
odep_bkp = opts.odep;

% generate nsmp samples by adding random noise to ss values
nsmp = 5;
[noisy_xss,noisy_fss] = addnoise(repmat(xss,1,nsmp),repmat(fss,1,nsmp));
[noisy_xdyn,noisy_fdyn] = addnoise_dyn(xdyn,fdyn,nsmp);


% perturb system from noisy initial conditions
pt_sol_id = [1 2 3 4 5 6];
[exp_sol,noisy_sol] =...
dopert_noisy_dyn(opts,noisy_xss,noisy_fss,odep_bkp,pt_sol_id,...
                noisy_xdyn,noisy_fdyn);
% [exp_sol,noisy_sol] = dopert_noisy(opts,noisy_xss,noisy_fss,odep_bkp,pt_sol_id);
close all