function [FX,DFX,D2FX,DFp,fx_sym,x,p] = obj_flux5_CAS
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
vemax = casadi.SX.sym('vemax',1);
KeFDP = casadi.SX.sym('KeFDP',1);
ne = casadi.SX.sym('ne',1);

% given fluxes
gflux = casadi.SX.sym('gflux',1);

flux = vemax.*(1-1./(1+(KeFDP./fdp).^ne)); % conv kin

x = [vemax;KeFDP;ne];
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
