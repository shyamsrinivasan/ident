function [FX,DFX,D2FX,DFp,fx_sym,x,p] = obj_flux2_k_CAS
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
K2pep = casadi.SX.sym('K2pep',1);
V2max = casadi.SX.sym('V2max',1);

% given fluxes
gflux = casadi.SX.sym('gflux',1);

% flux = k1cat.*enz.*ac./(ac+K1ac); % original model
flux = V2max.*pep./(pep+K2pep); % conv kin

x = [K2pep;V2max];
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
