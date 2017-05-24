%% no noise deterministic model 
% kinetic modeling with noisy metabolomics data 
% test using kotte model 
% algorithm :
% min |vmodel - vexpt|
%  st Svmodel = 0;
%     vmodel = f(x,p);
%     p >= pmin;
%     p <= pmax;
% given : x(expt), vexpt

%% load original kotte model
load('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\Kotte2014\model\kotte_model.mat');
p = pvec;
p(9) = [];
ival = [M;pvec(9)];
clear pvec

%% solve deterministic (no noise) model using casadi(cvodes)
odep = p;
solver_opts = struct('abstol',1e-3,'reltol',1e-3);
opts = struct('tspan',tspan,'x0',ival,'solver_opts',solver_opts,'odep',odep);
[xdyn,fdyn,xss1,fss1] = solveODE_cas(@kotte_CAS,opts,@kotte_flux_CAS);
figure
subplot(211);
plot(tspan,xdyn);
ylabel('concentrations a.u.');
legend('pep','fdp','E');
subplot(212)
plot(tspan,fdyn);
ylabel('fluxes a.u.');
legend('J','E(FDP)','vFbP','vEX','vPEPout');