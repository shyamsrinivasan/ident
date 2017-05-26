function [x4dyn,f4dyn,xss4,fss4] = check_kin_kotte(opts)
% check simulations from optimized CK model for kotte network

% solver_opts = struct('abstol',1e-6,'reltol',1e-6);
% opts = struct('tspan',tspan,'x0',xss1,'solver_opts',solver_opts,'odep',odep);
[x4dyn,f4dyn,xss4,fss4] =...
solveODE_cas(@kotte_CAS,opts,@kotte_flux_noCAS);
figure
subplot(211);
plot(opts.tspan,x4dyn);
ylabel('concentrations a.u.');
legend('pep','fdp','E','acetate');
subplot(212)
plot(opts.tspan,f4dyn);
ylabel('fluxes a.u.');
legend('J','E(FDP)','vFbP','vEX','vPEPout');