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
tspan = 0:0.1:500;
run_nonoise

% perturb model from ss
% perturb enzyme 2 (flux(4),vEX)
odep(14) = .5; % 2;
perturb_nonoise

%% noisy stochastic model - uses ode23 to solve
% tspan = 0:0.1:500;
run_withnoise

% perturb noisy model from ss
% perturb enzyme 2 (flux(4),vEX)
odep.p(14) = .5; % 2;
perturb_noise


%% casadi test
% [FX,FX_flx,fx_sym] = kotte_CAS();

% generate model data from convinience kinetics
[FX,~,~,fx_sym,~,~,FX_flx] = kotte_conkin_CAS();
clear tspan odep solver_opts opts
tspan = 0:0.1:2000;
p(14) = .5; % 2;
p(15) = .02;
p(16) = .3;
p(17) = 0;
odep = p;
solver_opts = struct('abstol',1e-3,'reltol',1e-3);
opts = struct('tspan',tspan,'x0',xss1(1:3),'solver_opts',solver_opts,'odep',odep);
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

% optimize for parameters by minimizing L2 norm of fluxes






