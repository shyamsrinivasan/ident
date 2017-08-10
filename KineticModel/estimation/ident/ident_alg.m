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
[exp_sol,noisy_sol] = dopert_noisy(opts,noisy_xss,noisy_fss,odep_bkp,pt_sol_id);
close all

ident_data = struct('nvar',3);

% 1. start by solving optimal estimation problem for fixed parameter i
modelfh = @(p)kotte_output(p,ident_data);
modelsensfh = @(p)kotte_sensoutput(p,data);
funh = struct('modelfh',modelfh,'modelsensfh',modelsensfh);

objfh = @(p)ident_obj(p,data,funh);
gradobjfh = @(p)ident_gradobj(p,data,funh);

prob = struct('obj',objfh,'gradobj',gradobjfh,'lb',lb,'ub',ub);
% optsol = nlconsopt(prob,data);

[xopt] = nlsqopt(prob,data);
% 2. take inc/dec. steps in parameter i
% 3. re-optimize parameters for new parameter i
% 4. repeat steps until threshold is crossed