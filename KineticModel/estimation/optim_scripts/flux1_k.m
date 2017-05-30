% optimization of flux parameters in kotte model for a CK formulation
function [x_opt,opt_id,new_opt_p,fval] =...
        flux1_k(opts,xss2,fss2,plist,old_opt_p)
if nargin<5
    old_opt_p = [];
end    

% flux 1    
p_id = cellfun(@(x)strcmpi(plist,x),{'K1ac','k1cat'},'UniformOutput',false);
p_id = cellfun(@(x)find(x),p_id);
% p = opts.odep(p_id)';
p = [.1;.1];

f2 = fss2(1); % add steayd state experimental flux
optim_p = [xss2;f2]; % concentrations & fluxes (expt) are parameters
lb = [1e-6;1e-3];
ub = [20;2000];
[x_opt,fval,~,~,opts] = runoptim_flux(opts,@obj_flux1_k_CAS,lb,ub,p,optim_p,1);

% check flux using conkin rate law
if ~isempty(old_opt_p)
    opts.odep = old_opt_p;
% else
%     pconv = [.1;.3;0]; % extra parameters for CK 'K1pep','K2fdp','rhoA'
%     opts.odep = [opts.odep';pconv];
end
opt_id = p_id; % [p_id,14];
opts.odep(opt_id) = x_opt;
new_opt_p = opts.odep;
