% call ode solver to solve ode
function [yout,flux] = solve_ode(fh,opts,flxh)

if isfield(opts,'tspan')
    tspan = opts.tspan;
end
if isfield(opts,'solver_opts')
    solver_opts =  opts.solver_opts;
else
    solver_opts = [];
end
if isfield(opts,'x0');
    x0 = opts.x0;
end
if isfield(opts,'odep')
    odep = opts.odep;
else
    odep = [];
end

if ~isempty(odep)
    newfh = @(t,x)fh(t,x,odep);
else
    newfh = fh;
end
tstart = tic;
[~,yout] = ode45(newfh,tspan,x0,solver_opts);
fprintf('\nTime to solve ode :%4.3f\n',toc(tstart));
yout = yout';

flux = flxh([yout;repmat(odep.model.PM,1,size(yout,2))],odep);