function solveODE_cas(opts)

if isfield(opts,'tspan')
    tspan = opts.tspan;
end
if isfield(opts,'solver_opts')
    solver_opts =  opts.solver_opts;
end
if isfield(opts,'x0');
    x0 = opts.x0;
end
if isfield(opts,'p')
    p = opts.p;
end

[tout
