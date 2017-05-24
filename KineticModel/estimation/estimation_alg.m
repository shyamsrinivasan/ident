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
run_nonoise

% perturb model from ss
% perturb enzyme 2 (flux(4),vEX)
odep(13) = .5; % 2;
perturb_nonoise

%% test - generate model data from convinience kinetics
% [FX,~,~,fx_sym,~,~,FX_flx] = kotte_conkin_CAS();
clear tspan odep solver_opts opts
tspan = 0:0.1:2000;
p(14) = 1; % 2;V2max
p(15) = .1; % K1pep
p(16) = .3; % K2fdp
p(17) = 0; % rhoA - regulation
odep = p;
solver_opts = struct('abstol',1e-3,'reltol',1e-3);
opts = struct('tspan',tspan,'x0',ival(1:3),'solver_opts',solver_opts,'odep',odep);
[x3dyn,f3dyn,xss3,fss3] = solveODE_cas(@kotte_conkin_CAS,opts,@kotte_convkinflux_CAS);
figure
subplot(211);
plot(tspan,x3dyn);
ylabel('concentrations a.u.');
legend('pep','fdp','E');
subplot(212)
plot(tspan,f3dyn);
ylabel('fluxes a.u.');
legend('J','E(FDP)','vFbP','vEX','vPEPout');

% optimization problem based on minimization of norm of difference in
% fluxes to estimate parameters for one flux
% get expression for objective function and grad(obj) for fluxi from casadi 
% [FX,DFX,D2FX,DFp,fx_sym,x_sym,p_sym,FX_flux] = obj_CAS();
[FX,gradF,~,DFp,fx_sym] = obj_flux1_CAS();
p1 = [p([1;11])';.1];
f2 = fss2(1);
optim_opts = [xss2;f2]; % steady state fluxes (expt) are parameters
obj = @(x)FX(x,optim_opts);
grad = @(x)gradF(x,optim_opts);
% contr = @(x)estimation_constr(x,optim_opts);
lb = [1e-6;1e-3;1e-6];
ub = [20;2000;20];
x0_opt = p1;

% setup opti problem
opt_prob = opti('obj',obj,'grad',grad,'bounds',lb,ub);
[x_opt,fval,exitflag,info] = solve(opt_prob,x0_opt); 



%% noisy stochastic model - uses ode23 to solve
% tspan = 0:0.1:500;
% run_withnoise

% perturb noisy model from ss
% perturb enzyme 2 (flux(4),vEX)
% odep.p(14) = .5; % 2;
% perturb_noise


% optimize for parameters by minimizing L2 norm of fluxes






