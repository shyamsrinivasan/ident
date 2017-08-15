% function to solve nlsq problem using casadi built in ipopt solver
function optsol = solve_nlsqopt(prob,x0)
if isfield(prob,'obj')
    obj = prob.obj;
end
if isfield(prob,'p')
    p = prob.p;
end    
if isfield(prob,'lb')
    lb = prob.lb;
else
    lb = zeros(length(x0),1);
end
if isfield(prob,'ub')
    ub = prob.ub;
else
    ub = zeros(length(x0),1);
    ub(ub==0) = Inf;
end

% setup nlp
nlp = struct('x',p,'f',obj);
solver = casadi.nlpsol('solver','ipopt',nlp);
% p0 = opts.odep(2:13)';
% p0(6) = .1;
sol_opt = solver('x0',x0,'lbx',lb,'ubx',ub);

optsol.xval = full(sol_opt.x);
optsol.fval = full(sol_opt.f);
optsol.lambda = full(sol_opt.lam_x);