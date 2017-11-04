function [FX,DFX,D2FX,DFp,fx_sym,x,p,FX_flux] = obj_CAS
% get function expression for objective function of optimization problem
% and its gradient

% x - parameters that are being optimized for
% p - parameters and concentrations that are being held constant

% concentrations
pep = casadi.SX.sym('pep',1);
fdp = casadi.SX.sym('fdp',1);
enz = casadi.SX.sym('enz',1);

% parameters
d = casadi.SX.sym('d',1);
ac = casadi.SX.sym('ac',1);
K1ac = casadi.SX.sym('K1ac',1);
K3fdp = casadi.SX.sym('K3fdp',1);
L3fdp = casadi.SX.sym('L3fdp',1);
K3pep = casadi.SX.sym('K3pep',1);
K2pep = casadi.SX.sym('K2pep',1);
vemax = casadi.SX.sym('vemax',1);
KeFDP = casadi.SX.sym('KeFDP',1);
ne = casadi.SX.sym('ne',1);
V4max = casadi.SX.sym('V4max',1);
k1cat = casadi.SX.sym('k1cat',1);
V3max = casadi.SX.sym('V3max',1);
V2max = casadi.SX.sym('V2max',1);

% given fluxes
gflux = casadi.SX.sym('gflux',5,1);

% change x and p such that constant concentrations are parameters and
% variable parameters are in x
x = [K1ac;K3fdp;L3fdp;K3pep;K2pep;vemax;KeFDP;ne;d;...
    V4max;k1cat;V3max;V2max];
p = [pep;fdp;enz;ac];

% fluxes - original model (for testing only) - subs w/ CK model
[flux,FX_flux] = kotte_flux_CAS(p,x);

diff_flux = [flux{1}-gflux(1);flux{2}-gflux(2);flux{3}-gflux(3);...
            flux{4}-gflux(4);flux{5}-gflux(5)];
p = [p;gflux];        
        
% objective function : min |v-v*|
fx_sym = norm(diff_flux);
FX = casadi.Function('FX',{x,p},{fx_sym});

% gradient of objectiv efunction
% get gradient of fluxes w.r.t parameters in x
% dflux_p = jacobian([flux{1};flux{2};flux{3};flux{4};flux{5}],x);
% grad = 1/(2*sqrt(fx_sym))*(sum(2.*diff_flux.*dflux_p));
% gradF = casadi.Function('gradF',{x,p},{grad});

dfx = jacobian(fx_sym,x);
DFX = casadi.Function('DFX',{x,p},{dfx});

dfp = jacobian(fx_sym,p);
DFp = casadi.Function('DFp',{x,p},{dfp});

D2FX = [];

% f(x) = 0
% fx_sym = cell(3,1);
% % PEP
% fx_sym{1} = flux1 - flux4 - flux5;
% % FBP
% fx_sym{2} = flux4 - flux3;
% % enzymes
% % E
% fx_sym{3} = flux2 - d*enz;

