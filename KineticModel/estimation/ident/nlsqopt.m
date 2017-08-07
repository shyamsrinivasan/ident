% optimization step for PLE-based identifiability
% idx - index of parameter that is fixed
function nlsqopt(x,p,data,idx)

% using casadi to solve nonlinear (objective) optimization 
% objective : y* - y
% bounds : p
% nlp = struct('x',var,'f',obj,'g',cons);
% solver = nlpsol('sol','ipopt',nlp);