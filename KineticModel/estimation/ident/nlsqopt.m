% optimization step for PLE-based identifiability
% idx - index of parameter that is fixed
function nlsqopt(x,p,prob,data)
if isfield(prob,'modelh')
    modelh = prob.modelh;
end
if isfield(prob,'cons')
    cons = prob.cons;
else
    cons = [];
end
if isfield(prob,'lb')
    lb = prob.lb;
end
if isfield(prob,'ub')
    ub = prob.ub;
end

% info
if isfield(data,'xexp')
    yexp = data.xexp;
end
if isifield(data,'xexp_var')
    yexp_var = data.yexp_var;
end
% create nlsqopt objective from modelh (modelh is a casadi fun 
% that gives a time course profile of y = f(y,p))
ymodel = modelh(x,p);
obj = sum((yexp-ymodel)./yexp_var).^2;


% using casadi to solve nonlinear (objective) optimization 
% objective : y* - y
% bounds : p
% nlp = struct('x',var,'f',obj,'g',cons);
% solver = nlpsol('sol','ipopt',nlp);