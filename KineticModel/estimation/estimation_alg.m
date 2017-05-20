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

odep = struct('p',p,'model',model);

% systems check
noisy_model = @(t,x)simnoisyODE_kotte(t,x,odep);
noisy_flux_test = noisyflux_kotte([M;model.PM],odep);

tspan = 0:0.1:200;
solver_opts = odeset('RelTol',1e-3,'AbsTol',1e-3);
% ac = find(strcmpi(model.mets,'ac[e]'));
% npts = 1;
% ap = 9;

% simulate noisy system from random initial condition
opts = struct('tspan',tspan,'x0',ival,'solver_opts',solver_opts,'odep',odep);
[xdyn,fdyn] = solve_ode(@simnoisyODE_kotte,opts,@flux_kotte);
figure
subplot(211);
plot(tspan,xdyn);
ylabel('concentrations a.u.');
subplot(212)
plot(tspan,fdyn);
ylabel('fluxes a.u.');

% perturb system from steady state and generate noisy data

% generate model data from convinience kinetics

% optimize for parameters by minimizing L2 norm of fluxes






