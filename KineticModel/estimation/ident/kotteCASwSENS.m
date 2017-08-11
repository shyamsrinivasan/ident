function [FXaug,FX,D2FX,fx_aug,aug_x,p] = kotteCASwSENS(nvar,np)

x = casadi.SX.sym('x',nvar,1);
p = casadi.SX.sym('p',np,1);
sens_var = casadi.SX.sym('sens_var',nvar,np);

fx_sys = cell(3,1);
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
fx_sens = [fx_sens{:}];
fx_sens_vec = reshape(fx_sens,[nvar*np,1]);
sens_var_vec = reshape(sens_var,[nvar*np,1]);

aug_x = [x;sens_var_vec];
fx_aug = [[fx_sys{:}]';fx_sens_vec(:)];

FXaug = casadi.Function('FXaug',{aug_x,p},{fx_aug});
FX = casadi.Function('FX',{x,p},{[fx_sys{:}]'});

% DFX = [];
D2FX = [];



% DFP = casadi.Function('DFP',{x,p},{dfp});


