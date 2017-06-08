% optimization of flux parameters in kotte model for a CK formulation
function [x_opt,opt_id,new_opt_p,fval] =...
        flux4_k(opts,xss2,fss2,plist,old_opt_p)
if nargin<5
    old_opt_p = [];
end    

% flux 1    
p_id = cellfun(@(x)strcmpi(plist,x),{'V4max'},'UniformOutput',false);
p_id = cellfun(@(x)find(x),p_id);
p = opts.odep(p_id)';% [.1;.5;.1;3e6]; 

f4 = fss2(5,:); % add steayd state experimental flux
optim_p = [xss2;...
           repmat(opts.odep(17),1,size(xss2,2));...
           f4]; % concentrations & fluxes (expt) are parameters
lb = 1e-3;
ub = 10;
[x_opt,fval,~,~,opts] = runoptim_flux(opts,@obj_flux4_k_CAS,lb,ub,p,optim_p,1);

% check flux using conkin rate law
if ~isempty(old_opt_p)
    opts.odep = old_opt_p;
end
opt_id = p_id; 
opts.odep(opt_id) = x_opt;
new_opt_p = opts.odep;
