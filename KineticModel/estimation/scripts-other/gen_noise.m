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
% dopert_noisy
% close all

% get only data from one steady state
pss = ones(length(pt_val.exp_pval),1);
pss(noisy_sol.xss(1,:)>noisy_sol.xss(2,:)) = 0;    

% use noisy data as perturbed data analoges to estimate parameters
% optimdata = struct('xss',{noisy_xss},...
%                    'fss',{noisy_fss});

% problem defn
optimdata = struct('nvar',5,'nc',3,'vexp',noisy_sol.fss(:,logical(pss)),...
                    'p_id',[1 11],'flxid',1);

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

prob = struct('obj',obj,'lb',lb,'ub',ub);
solveropt = struct('solver','scip');
optsol = nlconsopt(prob,x0,solveropt);




