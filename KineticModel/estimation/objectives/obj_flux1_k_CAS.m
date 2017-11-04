function [FX,DFX,D2FX,DFp,fx_sym,x,p] = obj_flux1_k_CAS(np)
% get function expression for objective function of optimization problem
% and its gradient

% x - parameters that are being optimized for
% p - parameters and concentrations that are being held constant

% concentrations
pep = casadi.SX.sym('pep',1,np);
fdp = casadi.SX.sym('fdp',1,np);
enz = casadi.SX.sym('enz',1,np);
ac = casadi.SX.sym('ac',1,np);

% parameters
K1ac = casadi.SX.sym('K1ac',1);
k1cat = casadi.SX.sym('k1cat',1);
% K1pep = casadi.SX.sym('K1pep',1);

% given fluxes
gflux = casadi.SX.sym('gflux',1,np);

% flux = k1cat.*enz.*ac./(ac+K1ac); % original model
flux = k1cat.*enz.*(ac./K1ac)./(1 + ac./K1ac); % conv kin

x = [K1ac;k1cat];
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
