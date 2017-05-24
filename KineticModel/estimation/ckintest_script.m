% convinience kinetic test script
% test for kotte network
% load original kotte model
load('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\Kotte2014\model\kotte_model.mat');
p = pvec;
ival = M;
clear pvec

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

