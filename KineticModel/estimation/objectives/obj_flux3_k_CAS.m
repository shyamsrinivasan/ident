function [FX,DFX,D2FX,DFp,fx_sym,x,p] = obj_flux3_k_CAS(np)
% get function expression for objective function of optimization problem
% and its gradient

% np - number of perturbations 
% x - parameters that are being optimized for
% p - parameters and concentrations that are being held constant

% concentrations
pep = casadi.SX.sym('pep',1,np);
fdp = casadi.SX.sym('fdp',1,np);
enz = casadi.SX.sym('enz',1,np);
ac = casadi.SX.sym('ac',1,np);

% parameters
K3pep = casadi.SX.sym('K3pep',1);
V3max = casadi.SX.sym('V3max',1);
K3fdp = casadi.SX.sym('K3fdp',1);
L3fdp = 4e6; % casadi.SX.sym('L3fdp',1);

% given fluxes
gflux = casadi.SX.sym('gflux',1,np);

ratio = 1+fdp./K3fdp;
flux = V3max.*(ratio-1).*(ratio).^3./...
            (ratio.^4+L3fdp.*(1+pep./K3pep).^(-4));

x = [K3pep;V3max;K3fdp]; % L3fdp];
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
