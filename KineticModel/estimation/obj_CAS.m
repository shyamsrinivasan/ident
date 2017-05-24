function [FX,DFX,D2FX,fx_sym,x,p,FX_flux] = obj_CAS()
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

% fluxes
flux1 = k1cat.*enz.*ac./(ac+K1ac);
flux2 = vemax.*(1-1./(1+(KeFDP./fdp).^ne));
ratio = 1+fdp./K3fdp;
flux3 = V3max.*(ratio-1).*(ratio).^3./...
            (ratio.^4+L3fdp.*(1+pep./K3pep).^(-4));
flux4 = V2max.*pep./(pep+K2pep);
flux5 = V4max*pep;

% f(x) = 0
fx_sym = cell(3,1);
% PEP
fx_sym{1} = flux1 - flux4 - flux5;
% FBP
fx_sym{2} = flux4 - flux3;
% enzymes
% E
fx_sym{3} = flux2 - d*enz;

% change x and p such that constant concentrations are parameters and
% variable parameters are in x
x = [K2pep;V2max];
p = [pep;fdp;enz;K1ac;K3fdp;L3fdp;K3pep;vemax;KeFDP;ne;ac;d;...
    V4max;k1cat;V3max];

FX = casadi.Function('FX',{x,p},{[fx_sym{1};fx_sym{2};fx_sym{3}]});
FX_flux = casadi.Function('FX_flux',{[flux1;flux2;flux3;flux4;flux5]});

dfx = jacobian([fx_sym{1};fx_sym{2};fx_sym{3}],x);
DFX = casadi.Function('DFX',{x,p},{dfx});

D2FX = [];
