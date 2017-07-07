% script to generate noisy data by adding noise to estimated steady state
% fluxes and concentrations in kotte model

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
opts.x0 = noisy_xss(:,1);
opts.tspan = 0:.1:300;
opts.odep = odep_bkp;
pt_val = struct('exp_pid',{11},...
                'exp_pval',{[.5;1.0;1.5;2]}); 
sol = getperturbations(pt_val,@perturb_nonoise,opts);

% add noise to perturbed data
noisy_sol = sol;
pt_p = repmat(odep_bkp,length(pt_val.exp_pval),1);
pt_p(:,pt_val.exp_pid) = pt_val.exp_pval;
[noisy_sol.xss,noisy_sol.fss] = addnoise(sol.xss,pt_p);
noisy_sol.xss = sol.xss + random(pd,3,size(sol.xss,2));

% get only data from one steady state
pss = ones(length(pt_val.exp_pval),1);
pss(noisy_sol.xss(1,:)>noisy_sol.xss(2,:)) = 2;    

% use noisy data as perturbed data analoges to estimate parameters
% optimdata = struct('xss',{noisy_xss},...
%                    'fss',{noisy_fss});

% problem defn
optimdata = struct('nvar',5,'nc',3,'vexp',noisy_fss,'p_id',[1 11]);

% set objective
obj = @(x)objnoisy(x,odep_bkp,optimdata);
% set constraints - no constraints
optimdata.eps = .5; % max permissible deviation from experimental value
% set bounds - concentration
lb = zeros(optimdata.nvar,1);
lb(1:optimdata.nc) = xss1.*(1-optimdata.eps);
ub = zeros(optimdata.nvar,1);
ub(1:optimdata.nc) = xss1.*(1+optimdata.eps);
% set bounds - parameter
lb(optimdata.nc+1:end) = [1e-3;1e-3];
ub(optimdata.nc+1:end) = [10;10];

% initial values for consrained nl(or quadratic?) optimization
x0 = [xss1;odep_bkp(optimdata.p_id)'];

optimopts = optiset('solver','scip',...
                      'maxiter',5000,...
                      'maxfeval',500000,...
                      'tolrfun',1e-6,...
                      'tolafun',1e-6,...
                      'display','final');   
optimopts = optiset(optimopts,'maxnodes',10000000);   
optim_prob = opti('obj',obj,'bounds',lb,ub,'options',optimopts);
[xval,fval,exitflag,info] = solve(optim_prob,x0); 





