% objective for identifiability analysis
% sum of square minimization of variance weighted error
function fx = ident_obj(x,p,data)

x = casadi.SX.sym('x',3,1);
p = casadi.SX.sym('p',17,1);

flux = kotte_flux_CAS(x,p);

fx_sym = cell(3,1);
fx_sym{1} = flux{1} - flux{4} - flux{5};
fx_sym{2} = flux{4} - flux{3};
fx_sym{3} = flux{2} - flux{6};

fx = casadi.Function('FX',{x,p},{[fx_sym{1};...
                                  fx_sym{2};...
                                  fx_sym{3}]});
                              






