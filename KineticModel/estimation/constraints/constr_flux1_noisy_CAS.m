function [FX_nlcon,DFX_nlcon,fx_nlcon] = constr_flux1_noisy_CAS
% objective function for estimation of kientic paramaters from noisy
% perturbation data
% concentrations are also variables

% np - number of perturbation whose data is used
fx_nlcon = cell(5,1);

% optimization variables
pep = casadi.SX.sym('pep',1);
fdp = casadi.SX.sym('fdp',1);
enz = casadi.SX.sym('enz',1);
ac = casadi.SX.sym('ac',1);
K1ac = casadi.SX.sym('K1ac',1);
k1cat = casadi.SX.sym('k1cat',1);
e = casadi.SX.sym('e',1);
xf = [pep;fdp;enz]; 

% optimization parameters
xss = casadi.SX.sym('xss',3,1);
fss = casadi.SX.sym('fss',1);

% flux parameters
K2pep = casadi.SX.sym('K2pep',1);
V2max = casadi.SX.sym('V2max',1);
K3pep = casadi.SX.sym('K3pep',1);
V3max = casadi.SX.sym('V3max',1);
K3fdp = casadi.SX.sym('K3fdp',1);
L3fdp = casadi.SX.sym('L3fdp',1);
V4max = casadi.SX.sym('V4max',1);
vemax = casadi.SX.sym('vemax',1);
KeFDP = casadi.SX.sym('KeFDP',1);
K1pep = casadi.SX.sym('K1pep',1);
K2fdp = casadi.SX.sym('K2fdp',1);
rhoA = casadi.SX.sym('rhoA',1);
ne = casadi.SX.sym('ne',1);
d = casadi.SX.sym('d',1);
pf = [K1ac;K3fdp;L3fdp;K3pep;K2pep;vemax;KeFDP;ne;d;...
        V4max;k1cat;V3max;V2max;K1pep;K2fdp;rhoA;ac];

% nl constraints
% constraint 1 : norm of flux for which parameters are estimated (1x1)
model_flux = kotte_flux_CAS(xf,pf);
model_flux{1} = model_flux{1} + e;    
fx_nlcon{1} = sqrt(sum((model_flux{1}-fss).^2));    

% constraint 2 : norm of all ss concentrations (1x1)
fx_nlcon{2} = sqrt(sum((xf-xss).^2));

% constraint 3 : ss for concentrations (mx1)
fx_nlcon(3:5) = kotte_ode(xf,pf,model_flux);

x = [xf;K1ac;k1cat;e];
p = [K3fdp;L3fdp;K3pep;K2pep;vemax;KeFDP;ne;d;...
     V4max;V3max;V2max;K1pep;K2fdp;rhoA;ac;xss;fss];

% all constraints : (m+2)x1
FX_nlcon = casadi.Function('FX_nlcon',{x,p},{[fx_nlcon{1},...
                                              fx_nlcon{2},...
                                              fx_nlcon{3:5}]'});
Dfx_sym = jacobian([fx_nlcon{1},...
                     fx_nlcon{2},...
                     fx_nlcon{3:5}]',x);
DFX_nlcon = casadi.Function('DFX_nlcon',{x,p},{Dfx_sym});


function dx = kotte_ode(x,p,flux)
dx = cell(3,1);

dx{1} = flux{1}-flux{4}-flux{5};
dx{2} = flux{4}-flux{3};
dx{3} = flux{2}-flux{6};
% dx{4} = x(4)-x(4);



