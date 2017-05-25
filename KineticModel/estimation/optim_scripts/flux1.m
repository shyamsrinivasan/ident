% optimization of flux parameters in kotte model for a CK formulation
function [x_opt,opt_pid,fss4,f4dyn,xss4,x4dyn,new_opt_p,fval] =...
        flux1(opts,xss2,fss2,plist,old_opt_p)
if nargin<5
    old_opt_p = [];
end    

% flux 1    
p_id = cellfun(@(x)strcmpi(plist,x),{'K1ac','k1cat'},'UniformOutput',false);
p_id = cellfun(@(x)find(x),p_id);
p = opts.odep(p_id)';

p = [p;.1]; % add K1pep to list of parameters
f2 = fss2(1); % add steayd state experimental flux
optim_p = [xss2;f2]; % concentrations & fluxes (expt) are parameters
lb = [1e-6;1e-3;1e-6];
ub = [20;2000;20];
[x_opt,fval,~,~,opts] = runoptim_flux(opts,@obj_flux1_CAS,lb,ub,p,optim_p);

% check flux using conkin rate law
if ~isempty(old_opt_p)
    opts.odep = old_opt_p;
else
    pconv = [.1;.3;0]; % extra parameters for CK 'K1pep','K2fdp','rhoA'
    opts.odep = [opts.odep';pconv];
end
opt_pid = [p_id,14];
opts.odep(opt_pid) = x_opt;
new_opt_p = opts.odep;
% solver_opts = struct('abstol',1e-6,'reltol',1e-6);
% opts = struct('tspan',tspan,'x0',xss1,'solver_opts',solver_opts,'odep',odep);
[x4dyn,f4dyn,xss4,fss4] = solveODE_cas(@kotte_conkin_CAS,opts,@kotte_convkinflux_noCAS);
figure
subplot(211);
plot(opts.tspan,x4dyn);
ylabel('concentrations a.u.');
legend('pep','fdp','E');
subplot(212)
plot(opts.tspan,f4dyn);
ylabel('fluxes a.u.');
legend('J','E(FDP)','vFbP','vEX','vPEPout');
