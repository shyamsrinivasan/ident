% optimization of flux parameters in kotte model for a CK formulation
function [x_opt,opt_id,new_opt_p,fval] =...
        flux1_k_noisy(opts,xss,fss,plist)
    
if isfield(opts,'opt_x0')
    x0 = opts.opt_x0;
end

% flux 1    
p_id = cellfun(@(x)strcmpi(plist,x),{'K1ac','k1cat'},'UniformOutput',false);
p_id = cellfun(@(x)find(x),p_id);

optim_p = opts.odep;
opts.p_id = p_id;
ss_val = [xss(:,1:4);...
          fss(1,1:4)];

% all concentrations as well as kientic parameters are variables
% lb = [pep;fdp;e;ac] - acetate is a equality constraint (fixed parameter)
lb = [0;0;0;1e-3;1e-3;0]; 
ub = [20;20;20;10;10;20];
% [x_opt,fval,~,~,opts] =...
% nlconstoptim_flux(opts,@obj_flux1_noisy_CAS,lb,ub,x0,...
%                   optim_p,ss_val,0,@constr_flux1_noisy_CAS); % linear objective
[x_opt,fval,~,~,opts] =...
nlconstoptim_flux_noCAS(opts,@obj_flux1_noisy,lb,ub,x0,...
                  optim_p,ss_val,0,@constr_flux1_noisy); % linear objective              

opt_id = p_id; % [p_id,14];
if ~isempty(x_opt)
    opts.odep(opt_id) = x_opt(4:end-1);
else
    opts.odep(opt_id) = NaN;
end
new_opt_p = opts.odep;
