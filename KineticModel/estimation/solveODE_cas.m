% ode solver wrapper for problems defined using casadi
function [yout,fout,yss,fss] = solveODE_cas(fh,opts,flxh)
if nargin>2
    givenflux = 1;
else
    givenflux = 0;    
end


if isfield(opts,'tspan')
    tspan = opts.tspan;
end
if isfield(opts,'solver_opts')
    solver_opts =  opts.solver_opts;
end
% if isfield(solver_opts,'abstol')
%     abstol = solver_opts.abstol;
% end
% if isfield(solver_opts,'reltol')
%     reltol = solver_opts.reltol;
% end
if isfield(opts,'x0')
    x0 = opts.x0;
end
if isfield(opts,'odep')
    odep = opts.odep;
end

% get expressions for function(s) and its jacobian
[~,~,~,fx,x,p] = fh();

% define casdi problem structure
neq = size(fx,1);
dx = casadi.SX(neq,1);
for i = 1:neq
    dx(i) = fx{i};
end
problem = struct('x',x,'ode',dx,'p',p);
solver_opts.output_t0 = 1;
solver_opts.print_stats = 1;
solver_opts.grid = tspan;
% options = struct('output_t0',1,'print_stats',1,'grid',tspan);
Fint = casadi.integrator('Fint','cvodes',problem,solver_opts);
sol = Fint('x0',x0,'p',odep);
yout = full(sol.xf);

if givenflux
    fout = flxh(yout,odep);
else
    fout = [];
end  

yss = yout(:,end);
if ~isempty(fout)
    fss = fout(:,end);
else
    fss = [];
end

% beta testing - foward sensitivity analysis using Casadi/SUNDIALS
% Fint_fwd = Fint.factory('Fint_fwd',{'x0','p','fwd:p'},{'fwd:xf'});
% res = Fint_fwd('x0',x0,'p',odep,'fwd_p',1);
% fwd_xf = full(res.fwd_xf);

