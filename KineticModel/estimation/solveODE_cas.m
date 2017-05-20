% ode solver wrapper for problems defined using casadi
function solveODE_cas(fh,opts)

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
    odep = opts.odep;
end

% get expressions for function(s) and its jacobian
[~,~,~,fx,x,p] = fh();

% define casdi problem structure
problem = struct('x',x,'ode',fx,'p',p);
options = struct('output_t0',1,'print_stats',1,'grid',tspan);
Fint = casadi.integrator('Fint','cvodes',problem,options);
sol = Fint('x0',x0,'p',odep);
yout = full(sol.xf);

