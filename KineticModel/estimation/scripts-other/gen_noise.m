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
pd = makedist('Uniform','lower',min(xss1),'upper',max(xss1));
met_noise = random(pd,3,nsmp);

pd = makedist('Uniform','lower',min(fss1),'upper',max(fss1));
flx_noise = random(pd,6,nsmp);

noisy_xss = repmat(xss1,1,nsmp)+met_noise;
noisy_fss = repmat(fss1,1,nsmp)+flx_noise;

% figure
% boxplot([noisy_fss';fss1']);
% figure
% boxplot([noisy_xss';xss1']);

% use noisy data as perturbed data analoges to estimate parameters
% optimdata = struct('xss',{noisy_xss},...
%                    'fss',{noisy_fss});

% problem defn
optimdata = struct('nvar',5,'nc',3,'vexp',noisy_fss,'p_id',[1 11]);

% set objective
obj = @(x)objnoisy(x,odep_bkp,optimdata);
% set constraints - no constraints
optimdata.eps = .1; % max permissible deviation from experimental value
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
                      'tolrfun',1e-8,...
                      'tolafun',1e-8,...
                      'display','final');   
optimopts = optiset(optimopts,'maxnodes',10000000);   
optim_prob = opti('obj',obj,'bounds',lb,ub,'options',optimopts);
[xval,fval,exitflag,info] = solve(optim_prob,x0); 





