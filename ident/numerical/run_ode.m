% simulate ode using either SUNDIALS or ode45
function solution = run_ode(odefun, solver_opts, fun_p, flux_fun)
if nargin<4
    flux_fun = [];
end
if nargin<3
    fun_p = [];
end

if ~isempty(fun_p)
    odefunh = @(t,x)odefun(t,x,fun_p);
else
    odefunh = odefun;
end
if isfield(solver_opts,'tspan')
    tspan = solver_opts.tspan;
else
    tspan = 0:.1:300;
end
x0 = solver_opts.x0;

% set ode45 options for rel and abs tol
if isfield(solver_opts,'reltol')
    reltol = solver_opts.reltol;
else
    reltol = 1e-6;
end
if isfield(solver_opts,'abstol')
    abstol = solver_opts.abstol;
else
    abstol = 1e-8;
end
options = odeset('RelTol',reltol,'AbsTol',abstol);

%call to ode45 to simulate ode
[tout, yout] = ode45(odefunh, tspan, x0, options);

solution.t = tout;
solution.y = yout;
solution.yss = yout(end,:)';
if ~isempty(flux_fun)
    solution.flux = flux_fun(yout', fun_p);
    solution.fss = solution.flux(end, :)';
end
return

