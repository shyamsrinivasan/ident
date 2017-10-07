% function to solve nlsq problem using casadi built in ipopt solver
function optsol = solve_nlsqopt(prob,x0)
if isfield(prob,'obj')
    obj = prob.obj;
end
if isfield(prob,'p')
    p = prob.p;
end    
if isfield(prob,'x')
    x = prob.x;
else 
    x = [];
end
if isfield(prob,'lb')
    if ~isempty(prob.lb)
        lb = prob.lb;
    else
        lb = zeros(length(x0),1);
        lb(lb==0) = .001;
    end
else
    lb = zeros(length(x0),1);
end
if isfield(prob,'ub') 
    if ~isempty(prob.ub)
        ub = prob.ub;
    else
        ub = zeros(length(x0),1);
        ub(ub==0) = Inf;
    end
else
    ub = zeros(length(x0),1);
    ub(ub==0) = Inf;
end

% setup nlp
if ~isempty(x)
    nlp = struct('x',p,'p',x,'f',obj);
else
    nlp = struct('x',p,'f',obj);
end
opts.ipopt.max_iter = 10000;
opts.ipopt.max_cpu_time = 1.5e6;
% opts.ipopt.fixed_variable_treatment = 'make_constraint';
solver = casadi.nlpsol('solver','ipopt',nlp,opts);
% p0 = opts.odep(2:13)';
% p0(6) = .1;
sol_opt = solver('x0',x0,'lbx',lb,'ubx',ub);

optsol.xval = full(sol_opt.x);
optsol.fval = full(sol_opt.f);
optsol.lambda = full(sol_opt.lam_x);
optsol.info = solver.stats;