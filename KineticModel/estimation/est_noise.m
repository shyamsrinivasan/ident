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

% get initial feasible solution for optimization based on these parameters 
% and fluxes

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
                   'fss',{sol(1).fss},...
                   'pid',{sol(1).exp_pid},...
                   'pval',{sol(1).exp_pval});

% set constraint rhs
m = 3; % concentrations
n = 1; % fluxes
nlrhs = zeros(2*m+n,1);
% flux norm
nlrhs(1:n) = 0.2;     
% concentration norm
nlrhs(n+1:n+m) = 0.3;
nlrhs(n+m+1:end) = 1e-8; %1e-6; % precision level for |Sv| <= 0
opts.nlrhs = nlrhs;

% constraint type : (-1 <=, 0 =, 1 >=)
nle = zeros(2*m+n,1);
nle(1:n) = -1;
nle(n+1:n+m) = -1;
% nle(n+m+1:n+2*m) = -1;
opts.nle = nle;

% opts.solver = 'nlopt';
% opts.opt_alg = 'GN_ISRES';
opts.solver = 'scip';
odep_opt = odep_bkp;

%% flux 1
opt_sol_f1 = runoptimp(opts,plist,odep_opt,optimopts,@optimize_flux1);   

%% flux 2
opt_sol_f2 = runoptimp(opts,plist,odep_opt,optimopts,@optimize_flux2);   

%% flux 3
opt_sol_f3 = runoptimp(opts,plist,odep_opt,optimopts,@optimize_flux3);  






% constrhs = constr_flux1_noisy(x0,opts.odep,[1,11],[sol.xss(:,5);sol.fss(:,5)]);








