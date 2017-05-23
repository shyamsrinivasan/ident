function [FX,DFX,D2FX,fx_sym,x,p,FX_flux] = kotte_conkin_CAS()

x = casadi.SX.sym('x',3,1);
p = casadi.SX.sym('p',17,1);

fx_sym = cell(3,1);

% parameters
d = p(10);

flux = kotte_convkinflux_CAS(x,p);
% PEP
fx_sym{1} = flux{1} - flux{4} - flux{5};
% FBP
fx_sym{2} = flux{4} - flux{3};
% enzymes
% E
fx_sym{3} = flux{2} - d*x(3);

FX = casadi.Function('FX',{x,p},{[fx_sym{1};...
                                  fx_sym{2};...
                                  fx_sym{3}]});
FX_flux = casadi.Function('FX_flux',{x,p},{[flux{1};...
                                            flux{2};...
                                            flux{3};...
                                            flux{4};...
                                            flux{5}]});
                                        
DFX = [];
D2FX = [];
                                        