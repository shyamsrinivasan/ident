function [x4dyn,f4dyn,xss4,fss4] = check_conkin_kotte(opts)
% check simulations from optimized CK model for kotte network

[x4dyn,f4dyn,xss4,fss4] =...
solveODE_cas(@kotte_conkin_CAS,opts,@kotte_convkinflux_noCAS);
figure
subplot(211);
plot(opts.tspan,x4dyn);
ylabel('concentrations a.u.');
legend('pep','fdp','E');
subplot(212)
plot(opts.tspan,f4dyn);
ylabel('fluxes a.u.');
legend('J','E(FDP)','vFbP','vEX','vPEPout');