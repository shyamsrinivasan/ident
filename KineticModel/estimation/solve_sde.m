function [yout,flux,yss,fss] = solve_sde(fh,gh,opts,flxh)
if nargin<3
    flxh = [];
end

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
[yout,w] = sde_euler(newfh,gh,tspan,x0,solver_opts);
fprintf('\nTime to solve ode :%4.3f\n',toc(tstart));
yout = yout';

yss = yout(:,end);
if ~isempty(flxh)
    flux = flxh(yout,odep);
else
    flux = [];
end

if ~isempty(flux)
    fss = flux(:,end);
else
    fss = [];
end