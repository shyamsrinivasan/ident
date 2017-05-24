function [x3dyn,f3dyn,xss3,fss3] = ckintest_script(tspan)
% convinience kinetic test script
% test for kotte network
% load original kotte model
load('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\Kotte2014\model\kotte_model.mat');
p = pvec;
ival = [M;p(9)]; % add acetate
p(9) = []; % remove acetate
clear pvec

% [FX,~,~,fx_sym,~,~,FX_flx] = kotte_conkin_CAS();
clear tspan odep solver_opts opts
p(13) = 1; % 2;V2max
p(14) = .1; % K1pep
p(15) = .3; % K2fdp
p(16) = 0; % rhoA - regulation
odep = p;
solver_opts = struct('abstol',1e-3,'reltol',1e-3);
opts = struct('tspan',tspan,'x0',ival,'solver_opts',solver_opts,'odep',odep);
[x3dyn,f3dyn,xss3,fss3] = solveODE_cas(@kotte_conkin_CAS,opts,@kotte_convkinflux_noCAS);
figure
subplot(211);
plot(tspan,x3dyn);
ylabel('concentrations a.u.');
legend('pep','fdp','E');
subplot(212)
plot(tspan,f3dyn);
ylabel('fluxes a.u.');
legend('J','E(FDP)','vFbP','vEX','vPEPout');

