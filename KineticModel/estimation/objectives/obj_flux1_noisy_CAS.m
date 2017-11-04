function [FXobj,DFXobj,fx_sym] = obj_flux1_noisy_CAS

x = casadi.SX.sym('x',6,1);
c = casadi.SX.sym('c',6,1);

fx_sym = c'*x;
FXobj = casadi.Function('FXobj',{x,c},{fx_sym});

dfx_sym = jacobian(fx_sym,x);
DFXobj = casadi.Function('DFXobj',{x,c},{dfx_sym});