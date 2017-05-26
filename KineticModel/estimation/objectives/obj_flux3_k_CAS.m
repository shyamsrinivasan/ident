function [FX,DFX,D2FX,DFp,fx_sym,x,p] = obj_flux3_k_CAS
% get function expression for objective function of optimization problem
% and its gradient

% x - parameters that are being optimized for
% p - parameters and concentrations that are being held constant

% concentrations
pep = casadi.SX.sym('pep',1);
fdp = casadi.SX.sym('fdp',1);
enz = casadi.SX.sym('enz',1);
ac = casadi.SX.sym('ac',1);

% parameters
K3pep = casadi.SX.sym('K3pep',1);
V3max = casadi.SX.sym('V3max',1);
K3fdp = casadi.SX.sym('K3fdp',1);
L3fdp = casadi.SX.sym('L3fdp',1);

% given fluxes
gflux = casadi.SX.sym('gflux',1);

ratio = 1+fdp./K3fdp;
flux{3} = V3max.*(ratio-1).*(ratio).^3./...
            (ratio.^4+L3fdp.*(1+pep./K3pep).^(-4));

x = [K3pep;V3max;K3fdp;L3fdp];
p = [pep;fdp;enz;ac;gflux];

fx_sym = norm(flux-gflux);
FX = casadi.Function('FX',{x,p},{fx_sym});

% gradient
dfx = jacobian(fx_sym,x);
DFX = casadi.Function('DFX',{x,p},{dfx});

% parameter sensitivity
dfp = jacobian(fx_sym,p);
DFp = casadi.Function('DFp',{x,p},{dfp});

D2FX = [];
