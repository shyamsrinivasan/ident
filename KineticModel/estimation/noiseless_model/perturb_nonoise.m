function [x2dyn,f2dyn,xss2,fss2] = perturb_nonoise(opts)
% perturb model with no noise
% perturb system from steady state and generate data w/o noise
% perturb enzyme 2 (flux(4),vEX)
% opts.x0 = xss1(1:4);
% opts.odep = odep;
% [x2dyn,f2dyn,xss2,fss2] = solve_ode(@simnoisyODE_kotte,opts,@flux_kotte);
[x2dyn,f2dyn,xss2,fss2] = solveODE_cas(@kotte_CAS,opts,@kotte_flux_noCAS);
figure
subplot(211);
plot(opts.tspan,x2dyn);
ylabel('concentrations a.u.');
legend('pep','fdp','E');
subplot(212)
plot(opts.tspan,f2dyn);
ylabel('fluxes a.u.');
legend('J','E(FDP)','vFbP','vEX','vPEPout');

