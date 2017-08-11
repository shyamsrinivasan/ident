% function to provide time course for x and dx/dp
% p - parameter vector
% data - other required info structure
function yout = kotte_senseoutput(p,data)
if isfield(data,'nc')
    nc = data.nc;
end
if isfield(data,'np')
    np = data.np;
end    
if isfield(data,'x0')
    x0 = data.x0;
end

fh = @()kotteCASwSENS(nc,np); % @kotte_CAS
ival = [x0;zeros(nc*np,1)];

solver_opts = struct('abstol',1e-6,'reltol',1e-6);
opts = struct('tspan',data.tspan,'x0',ival,'solver_opts',solver_opts,'odep',p);
yout = solveODE_cas(fh,opts);