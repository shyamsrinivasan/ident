% function to provide time course for given x and p inside nlsqopt
% x - states
% p - parameter vector
% data - other required info structure
function yout = kotte_output(y,p,data)

solver_opts = struct('abstol',1e-6,'reltol',1e-6);
opts = struct('tspan',data.tspan,'x0',y,'solver_opts',solver_opts,'odep',p);
yout = solveODE_cas(@kotte_CAS,opts,@kotte_flux_noCAS);


