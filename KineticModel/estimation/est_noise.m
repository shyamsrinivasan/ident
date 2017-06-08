% script for estimating parameters from noisy data
plist = {'K1ac','K3fdp','L3fdp','K3pep','K2pep','vemax','KeFDP','ne',...
        'd','V4max','k1cat','V3max','V2max','K1pep','K2fdp','rhoA'};

% generate experimental data - get initial ss
tspan = 0:0.1:300;
[xdyn,fdyn,xss1,fss1,opts] = run_withnoise(tspan);   

% backup parameters and initial conditions
ival_bkp = opts.x0;
odep_bkp = opts.odep;

%% perturbation to all fluxes 
opts.x0 = xss1;
opts.tspan = 0:.1:600;
opts.odep = odep_bkp;
% flux 1, 2 and 3 % k1cat, 'V2max', 'V3max'
ptopts = struct('exp_pid',{11},...
                'exp_pval',{[.1;.5;1.0;1.5;2]}); 
sol = getperturbations(ptopts,@perturb_noisy,opts);

% run perturbation without noise to get feasible initial solution
tspan = 0:0.1:300;
[~,~,xss1i,fss1i,optsi] = run_nonoise(tspan);   
optsi.odep = odep_bkp;
optsi.x0 = xss1i;
optsi.tspan = 0:.1:600;
soli = getperturbations(ptopts,@perturb_nonoise,optsi);
% close all

% flux 4 V4max
% opts.tspan = 0:.1:10000;
% ptopts = struct('exp_pid',10,'exp_pval',[0;.1;.3;.5;.7;.9;1]);
% sol = getperturbations(ptopts,@perturb_noisy,opts,sol);
% close all

%% use single perturbation sets to get parameters
optimopts = struct('xss',{sol(1).xss},...
                   'fss',{sol(1).fss});
opts.init_xss = soli(1).xss(:,1);   
opts.solver = 'nlopt';
% opts.opt_alg = 'GN_ISRES';
odep_opt = odep_bkp;
% odep_opt(11) = 2;
opt_sol = runoptimp(opts,plist,odep_opt,optimopts,@optimize_p_noisy);     





