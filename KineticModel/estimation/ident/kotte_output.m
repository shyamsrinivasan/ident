% function to provide time course for given x and p inside nlsqopt
% x - states
% p - parameter vector
% data - other required info structure
function yout = kotte_output(p,data)
if isfield(data,'x0')
    x0 = data.x0;
end

solver_opts = struct('abstol',1e-6,'reltol',1e-6);
opts = struct('tspan',data.tspan,'x0',x0,'solver_opts',solver_opts,'odep',p);
yout = solveODE_cas(@kotte_CAS,opts);