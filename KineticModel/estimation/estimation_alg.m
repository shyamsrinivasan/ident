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

% load original kotte model
load('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\Kotte2014\model\kotte_model.mat');
p = pvec;
ival = M;
clear pvec

% solve using casadi
tspan = 0:0.1:500;
% odep = p;
% solver_opts = struct('abstol',1e-3,'reltol',1e-3);
% opts = struct('tspan',tspan,'x0',ival,'solver_opts',solver_opts,'odep',odep);
% [xdyn,fdyn,xss1,fss1] = solveODE_cas(@kotte_CAS,opts,@kotte_flux_CAS);

odep = struct('p',p,'model',model);

% systems check
noisy_model = @(t,x)simnoisyODE_kotte(t,x,odep);
noisy_flux_test = noisyflux_kotte([M;model.PM],odep);

solver_opts = odeset('RelTol',1e-3,'AbsTol',1e-3);
opts = struct('tspan',tspan,'x0',ival,'solver_opts',solver_opts,'odep',odep);

% simulate noisy system from random initial condition - this is better for
% simulating the noisy model
[xdyn,fdyn,xss1,fss1] = solve_ode(@simnoisyODE_kotte,opts,@flux_kotte);

figure
subplot(211);
plot(tspan,xdyn);
ylabel('concentrations a.u.');
legend('pep','fdp','E');
subplot(212)
plot(tspan,fdyn);
ylabel('fluxes a.u.');
legend('J','E(FDP)','vFbP','vEX','vPEPout');

% perturb system from steady state and generate noisy data
% perturb enzyme 2 (flux(4),vEX)
odep.p(14) = .5; % 2;
odep.p(15) = .02;
odep.p(16) = .3;
opts.ival = xss1(1:3);
opts.odep = odep;
[x2dyn,f2dyn,xss2,fss2] = solve_ode(@simnoisyODE_kotte,opts,@flux_kotte);
figure
subplot(211);
plot(tspan,xdyn);
ylabel('concentrations a.u.');
legend('pep','fdp','E');
subplot(212)
plot(tspan,fdyn);
ylabel('fluxes a.u.');
legend('J','E(FDP)','vFbP','vEX','vPEPout');

% casadi test
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






