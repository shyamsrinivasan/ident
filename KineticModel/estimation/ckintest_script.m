% convinience kinetic test script
% test for kotte network
% load original kotte model
load('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\Kotte2014\model\kotte_model.mat');
p = pvec;
ival = M;
clear pvec

p(15) = .02;
p(16) = .3;
p(17) = 0;

odep = struct('p',p,'model',model);

% systems check
convkin_model = @(t,x)conkinODE_kotte(t,x,odep);
convkin_flux_test = convkin_kotte([M;model.PM],odep);

tspan = 0:0.1:60000;
solver_opts = odeset('RelTol',1e-3,'AbsTol',1e-3);

% simulate system from random initial condition
opts = struct('tspan',tspan,'x0',ival,'solver_opts',solver_opts,'odep',odep);
[xdyn,fdyn,xss1,fss1] = solve_ode(@conkinODE_kotte,opts,@convkin_kotte);
figure
subplot(211);
plot(tspan,xdyn);
ylabel('concentrations a.u.');
legend('pep','fdp','E');
subplot(212)
plot(tspan,fdyn);
ylabel('fluxes a.u.');
legend('J','E(FDP)','vFbP','vEX','vPEPout');

