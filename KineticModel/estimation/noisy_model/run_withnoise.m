%% noisy stochastic model
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

%% load original kotte model
load('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\Kotte2014\model\kotte_model.mat');
p = pvec;
ival = M;
clear pvec odep solver_opts opts

%% solve noisy model using ode45
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