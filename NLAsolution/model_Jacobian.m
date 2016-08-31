function J = model_Jacobian(model,mc,rhshandle,Jpatt)
if nargin < 4
    % get model sparsity pattern
    Jpatt = modelSparsity(model);
end

% obtain handle for for finite diference approximations
fDh = finiteD;

% obtain function handle for f(X) in Jacobian calculation
% fJrow = @(i,x)Svrow(i,x,model,pvec);

% execute finite difference approximation of jacobian with fh
J = fDh(Jpatt,mc,rhshandle);


% nested function
function fi = fcall(rowid,x)
% calculate rxnid
% rxnid <- rxns in which x(rowid) participates
fi = fJrow(model,pvec,x,rowid,rxnid);
% fi = feval(fi,[x1 x2..xn],rxnid)
% rxnid <- rxns in which xi participates
% flux = iflux(model,pvec,x,[rxnid])
% fi = S(i,:)*

function fi = fJrow(model,pvec,mc,rowid,rxnid)
flux = zeros(size(model.S,2),1);
flux(rxnid) = iflux(model,pvec,mc,flux,rxnid);

fi = model.S(rowid,:)*flux;