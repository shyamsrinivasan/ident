% solve ode for model for given parameter vector to get time course
% data using casadi built in cvodes
function yout = solve_model(aug_p,data)

if isfield(data,'ival')
    ival = data.ival;
end
if isifield(data,'tspan')
    tspan = data.tspan;
end
if isfield(data,'odefh')
    fh = data.odefh;
else
    fh = [];
end

% get casadi symbolics
if ~isempty(fh)
    [~,~,~,fx,xsym,psym] = fh();
    neq = size(fx,1);
    dx = casadi.SX(neq,1);
    for i = 1:neq
        dx(i) = fx{i};
    end
    problem = struct('x',xsym,'ode',dx,'p',psym);
    solver_opts.output_t0 = 1;
    solver_opts.print_stats = 1;
    solver_opts.grid = tspan;
    Fint = casadi.integrator('Fint','cvodes',problem,solver_opts);
    sol = Fint('x0',ival,'p',aug_p);
    yout = full(sol.xf);
end




