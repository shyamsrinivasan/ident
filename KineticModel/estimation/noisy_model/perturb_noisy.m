function [xdyn,fdyn,xss,fss] = perturb_noisy(opts)
% perturb model with no noise
% perturb system from steady state and generate data w/ noise
solver_opts = sdeset('SDEType','Ito',...
                 'RandSeed',2,...
                 'ConstGFUN','yes',...
                 'NonNegative','yes',...
                 'DiagonalNoise','yes'); % default for DiagonalNoise
opts.solver_opts = solver_opts;             
% g = eye(3);
g = 0.01*ones(3,1);

% [xdyn,fdyn,xss,fss] = solve_ode(@simnoisyODE_kotte,opts,@kotte_flux_noCAS);
[xdyn,fdyn,xss,fss] = solve_sde(@simnoisyODE_kotte,g,opts,@kotte_flux_noCAS);
figure
subplot(211);
plot(opts.tspan,xdyn);
ylabel('concentrations a.u.');
legend('pep','fdp','E');
subplot(212)
plot(opts.tspan,fdyn);
ylabel('fluxes a.u.');
legend('J','E(FDP)','vFbP','vEX','vPEPout','vEout');
