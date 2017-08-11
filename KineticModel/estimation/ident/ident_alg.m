% identifiability analysis algorithm
% generate noisy experimental data
% generate experimental data - get initial ss
tspan = 0:0.1:100;
[xdyn,fdyn,xss,fss,opts] = run_nonoise(tspan);

% backup parameters and initial conditions
ival_bkp = opts.x0;
odep_bkp = opts.odep;

% generate nsmp samples by adding random noise to ss values
nsmp = 10;
[noisy_xss,noisy_fss] = addnoise(repmat(xss,1,nsmp),odep_bkp);

% perturb system from noisy initial conditions
pt_sol_id = [1 2 3];
[exp_sol,noisy_sol] =...
dopert_noisy_dyn(opts,noisy_xss,noisy_fss,odep_bkp,pt_sol_id);
close all

xexp_var = ones(size(noisy_sol(1).xdyn,1),size(noisy_sol(1).xdyn,2));
optim_opts = struct('pname','K1ac','nc',3,'nf',6,...
                    'odep',odep_bkp,...
                    'tspan',0:.1:300,...
                    'x0',noisy_xss(:,1),...
                    'xexp_dyn',noisy_sol(1).xdyn,...
                    'xexp_var',xexp_var,...
                    'objf','ident_obj',...
                    'gradobjfh','ident_gradobj',...
                    'boundfh','boundsallp',...
                    'modelf','kotte_output',...
                    'modelsensf','kotte_senseoutput');

[prob,data] = setup_ident_prob(optim_opts);

% set initial parameter vales for optimization
p0 = opts.odep(data.varid)';
solveropt = struct('solver','ipopt','multi',0);

% 1. start by solving optimal estimation problem for fixed parameter i
[xopt] = nlsqopt(prob,p0,solveropt);
                
% 2. take inc/dec. steps in parameter i
% 3. re-optimize parameters for new parameter i
% 4. repeat steps until threshold is crossed