% estimate fluxes for flux 3 in kotte model fdp ---> 

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
[exp_sol,noisy_sol] = dopert_noisy(opts,noisy_xss,noisy_fss,odep_bkp,pt_sol_id);
close all

% get only data from one steady state
pss = ones(1,numel(exp_sol.exp_pval));
% pss(exp_sol.xss(1,:)>exp_sol.xss(2,:)) = 0;    

% options structure for solving problem
optimdata = struct('nc',3,'nflx',6,'nf',1,'flxid',3,'eps_v',.05,...
                    'eps_c',.3,'vexp',exp_sol.fss(:,logical(pss)),...
                    'xexp',exp_sol.xss(:,logical(pss)),...
                    'odep',odep_bkp,...
                    'wt_xss',xss(:,1),'wt_fss',fss(:,1));
                
% define all bound fun handles
allboundh = {'boundsf1',...
             [],...
             'boundsf3',...
             [],...
             [],...
             []};
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
            'consnoisyf3',...
            'consnoisyf2',...
            'consnoisyf1',...
            []};
% define all rhs for nl cons
allnlrhs = {0,[],[],0,0,[]};
% define cons type for nl cons
allnle = {0,[],[],0,0,[]};

% problem defn
setup_opts = optimdata;
setup_opts.bounds = allboundh;
setup_opts.obj = allobjh;
setup_opts.nlcons = allconsh;
setup_opts.nlrhs = allnlrhs;
setup_opts.nle = allnle;

[prob,optimdata] = setup_optim_prob(setup_opts);

% initial values for consrained nl(or quadratic?) optimization
x0 = getrandomivals(optimdata,.1,5000);
solveropt = struct('solver','ipopt','multi',0);
optsol = choose_nlconsopt(prob,x0,optimdata,solveropt);

% combine results for comparison plot
est_data = combine_results(optsol,opts,noisy_sol,optimdata,pss,pss);

% compare fluxes and concentrations
compare_vals(est_data,noisy_sol,optimdata,opts,pss);

% compare parameters in parameter space
compare_pars(est_data);




