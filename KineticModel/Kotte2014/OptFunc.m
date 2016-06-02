function out = OptFunc(x,model,pvec)
out{1} = @(x)objval(x,model,pvec);
out{2} = @(x)objgrad(x,model,pvec);
out{3} = @(x)nlconstr(x,model,pvec);
out{4} = @(x)nlconstrjac(x,model,pvec);
out{5} = @()nlconstrjacsparsity(model);

function z = objval(x,model,pvec)
mc = x(1:3);
pvec = cons([pvec x(4:end)],x);
if ~isempty(model)
    PM = cons(model.PM,mc);
    allmc = [mc;PM];
else
    allmc = mc;
end
flux = Kotte_givenFlux(allmc,pvec,model);
pepout = 5;
z = zeros(1,1);
z = cons(z,allmc);
mss_pepout = 1.9;
mss_pepout = cons(mss_pepout,allmc);
z(1) = flux(pepout) - mss_pepout;

function grad = objgrad(x,model)
xADobj = deriv(x,eye(length(x)));
ADfun = @(x)objval(x,model);
xADres = ADfun(xADobj);
grad = getydot(xADres);

function A = nlconstr(x,model)
mc = x(1:3);
pvec = x(4:end);
A = Kotte_givenNLAE(mc,model,pvec);

function jac = nlconstrjac(x,model)
xADobj = deriv(x,eye(length(x)));
ADfun = @(x)nlconstr(x,model);
xADres = ADfun(xADobj);
jac = getydot(xADres); 

function jspr = nlconstrsparsity(model)
Jpatt = modelSparsity(model);
Jmc = Jpatt(1:3,1:3);
% Jpvec = sparse([



