% perturb noisy model
% perturb system from steady state and generate noisy data
% perturb enzyme 2 (flux(4),vEX)
opts.x0 = xss1(1:3);
opts.odep = odep;
[x2dyn,f2dyn,xss2,fss2] = solve_ode(@simnoisyODE_kotte,opts,@flux_kotte);
figure
subplot(211);
plot(tspan,x2dyn);
ylabel('concentrations a.u.');
legend('pep','fdp','E');
subplot(212)
plot(tspan,f2dyn);
ylabel('fluxes a.u.');
legend('J','E(FDP)','vFbP','vEX','vPEPout');