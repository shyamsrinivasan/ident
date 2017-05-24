% kinetic modeling with noisy metabolomics data 
% test using kotte model 
% algorithm :
% min |vmodel - vexpt|
%  st Svmodel = 0;
%     vmodel = f(x,p);
%     p >= pmin;
%     p <= pmax;
% given : x(expt), vexpt

% experimental data is generated through addition of noise to actual kotte model
% vmodel will use convenience kinetics for estimating fluxes

%% no noise deterministic model - uses casadi to solve
tspan = 0:0.1:300;
[xdyn,fdyn,xss1,fss1,opts] = run_nonoise(tspan);

% perturb model from ss
% perturb enzyme 2 (flux(4),vEX)
opts.x0 = xss1;
opts.odep(13) = .5; % 2;
[x2dyn,f2dyn,xss2,fss2] = perturb_nonoise(opts);

%% test - generate model data from convinience kinetics
% tspan = 0:0.1:2000;
% [x3dyn,f3dyn,xss3,fss3] = ckintest_script(tspan);
%% run flux parameter estimation
% optimization problem based on minimization of norm of difference in
% fluxes to estimate parameters for one flux
% get expression for objective function and grad(obj) for fluxi from casadi 
% [FX,DFX,D2FX,DFp,fx_sym,x_sym,p_sym,FX_flux] = obj_CAS();
% list of all available parameters in order of occurence
plist = {'K1ac','K3fdp','L3fdp','K3pep','K2pep','vemax','KeFDP','ne',...
        'd','V4max','k1cat','V3max','V2max','K1pep','K2fdp','rhoA'};
p_id = cellfun(@(x)strcmpi(plist,x),{'K1ac','k1cat'},'UniformOutput',false);
p_id = cellfun(@(x)find(x),p_id);
p = opts.odep(p_id)';

p = [p;.1]; % add K1pep to list of parameters
f2 = fss2(1); % add steayd state experimental flux
optim_p = [xss2;f2]; % concentrations & fluxes (expt) are parameters
lb = [1e-6;1e-3;1e-6];
ub = [20;2000;20];
[x_opt,fval,~,~,opts] = runoptim_flux(opts,@obj_flux1_CAS,lb,ub,p,optim_p);

%% check flux using conkin rate law
pconv = [.1;.3;0]; % extra parameters for CK 'K1pep','K2fdp','rhoA'
odep = [opts.odep';pconv];
opt_pid = [p_id,14];
odep(opt_pid) = x_opt;
solver_opts = struct('abstol',1e-6,'reltol',1e-6);
opts = struct('tspan',tspan,'x0',xss1,'solver_opts',solver_opts,'odep',odep);
[x4dyn,f4dyn,xss4,fss4] = solveODE_cas(@kotte_conkin_CAS,opts,@kotte_convkinflux_noCAS);
figure
subplot(211);
plot(tspan,x4dyn);
ylabel('concentrations a.u.');
legend('pep','fdp','E');
subplot(212)
plot(tspan,f4dyn);
ylabel('fluxes a.u.');
legend('J','E(FDP)','vFbP','vEX','vPEPout');




%% noisy stochastic model - uses ode23 to solve
% tspan = 0:0.1:500;
% run_withnoise

% perturb noisy model from ss
% perturb enzyme 2 (flux(4),vEX)
% odep.p(14) = .5; % 2;
% perturb_noise


% optimize for parameters by minimizing L2 norm of fluxes






