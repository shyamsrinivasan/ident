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
                'exp_pval',{2}); 
sol = getperturbations(ptopts,@perturb_noisy,opts);
% close all

% flux 4 V4max
% opts.tspan = 0:.1:10000;
% ptopts = struct('exp_pid',10,'exp_pval',[0;.1;.3;.5;.7;.9;1]);
% sol = getperturbations(ptopts,@perturb_noisy,opts,sol);
% close all





