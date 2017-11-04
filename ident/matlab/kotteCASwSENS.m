function [FXaug,FX,D2FX,fx_aug,aug_x,p] = kotteCASwSENS(nvar,np,idx)

x = casadi.SX.sym('x',nvar,1);
p = casadi.SX.sym('p',np,1);
sens_var = casadi.SX.sym('sens_var',nvar,np);

fx_sys = cell(nvar,1);
fx_sens = cell(np,1);

% system equations
flux = kotte_flux_CAS(x,p);
fx_sys{1} = flux{1} - flux{4} - flux{5};
fx_sys{2} = flux{4} - flux{3};
fx_sys{3} = flux{2} - flux{6};

% system jacobian
dfx = jacobian([fx_sys{:}],x);

% parameter senstivity equations
dfp = jacobian([fx_sys{:}],p);

for ip = 1:np
    fx_sens{ip} = dfx*sens_var(:,ip) + dfp(:,ip);
end
if ~isempty(idx)
    pid = setdiff(1:np,idx);    
else
    pid = 1:np;    
end
fx_sens = [fx_sens{pid}]; % remove column corresponding to idx parameter
fx_sens_vec = reshape(fx_sens,[nvar*length(pid),1]);
sens_var_vec = reshape(sens_var(:,pid),[nvar*length(pid),1]);

aug_x = [x;sens_var_vec];
fx_aug = [[fx_sys{:}]';fx_sens_vec(:)];

FXaug = casadi.Function('FXaug',{aug_x,p},{fx_aug});
FX = casadi.Function('FX',{x,p},{[fx_sys{:}]'});

% DFX = [];
D2FX = [];



% DFP = casadi.Function('DFP',{x,p},{dfp});


