% script to estimate parameters from noisy data for flux 1 in kotte model

% generate experimental data - get initial ss
tspan = 0:0.1:100;
[xdyn,fdyn,xss1,fss1,opts] = run_nonoise(tspan);

% backup parameters and initial conditions
ival_bkp = opts.x0;
odep_bkp = opts.odep;

% generate nsmp samples by adding random noise to ss values
nsmp = 10;
[noisy_xss,noisy_fss] = addnoise(repmat(xss1,1,nsmp),odep_bkp);

% figure
% boxplot([noisy_fss';fss1']);
% figure
% boxplot([noisy_xss';xss1']);

% perturb system from noisy initial conditions
pt_sol_id = [1 2 3];
[exp_sol,noisy_sol] = dopert_noisy(opts,noisy_xss,odep_bkp,pt_sol_id);
close all

% get only data from one steady state
pss = ones(1,numel(exp_sol.exp_pval));
% pss(exp_sol.xss(1,:)>exp_sol.xss(2,:)) = 0;    

% options structure for solving problem
optimdata = struct('nc',3,'nf',1,'flxid',1,'eps',.5,...
                    'vexp',exp_sol.fss(:,logical(pss)),...                    
                    'odep',odep_bkp,...
                    'wt_xss',noisy_xss(:,1),'wt_fss',noisy_fss(:,1));

% define all objectiv efun handles
allobjh = {'objnoisy',...
           [],...
           'objnoisy',...
           'objnoisy',...
           'objnoisy',...
           []};
% define all constraint fun handles
allconsh = {'consnoisyf1',...
            [],...
            'consnoisyf1',...
            'consnoisyf1',...
            'consnoisyf1',...
            []};
% define all rhs for nl cons
allnlrhs = {0,[],0,0,0,[]};
% define cons type for nl cons
allnle = {0,[],0,0,0,[]};

% problem defn
setup_opts = optimdata;
setup_opts.obj = allobjh;
setup_opts.nlcons = allconsh;
setup_opts.nlrhs = allnlrhs;
setup_opts.nle = allnle;

[prob,optimdata] = setup_optim_prob(setup_opts);

% initial values for consrained nl(or quadratic?) optimization
x0 = [noisy_xss(:,2);optimdata.odep(optimdata.p_id)';noisy_fss(1,2)];

% prob = struct('obj',obj,'nlcons',nlcons,'nlrhs',nlrhs,'nle',nle,'lb',lb,'ub',ub);
solveropt = struct('solver','ipopt','multi',1,'multi_pts',[7 7]);
optsol = nlconsopt(prob,x0,solveropt,optimdata);

% compare fluxes and concentrations
compare_vals(optsol,noisy_sol,optimdata,opts);




