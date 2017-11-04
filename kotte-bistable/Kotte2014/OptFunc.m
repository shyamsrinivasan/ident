function out = OptFunc(model,pvec)
out{1} = @(x)objval(x,model,pvec);
out{2} = @(x)objgrad(x,model,pvec);
out{3} = @(x)nlconstr(x,model,pvec);
out{4} = @(x)nlconstrjac(x,model,pvec);
out{5} = @()nlconstrjacsparsity(model,pvec);

function z = objval(x,model,pvec)
mc = x(1:3);
pvec = cons([pvec x(4:end)'],x);
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
z(1) = (flux(pepout) - mss_pepout)^2;

function grad = objgrad(x,model,pvec)
xADobj = deriv(x,eye(length(x)));
ADfun = @(x)objval(x,model,pvec);
xADres = ADfun(xADobj);
grad = getydot(xADres);

function A = nlconstr(x,model,pvec)
mc = x(1:3);
pvec = cons([pvec x(4:end)'],x);
A = Kotte_givenNLAE(mc,model,pvec);

function jac = nlconstrjac(x,model,pvec)
% method 1 - jacobian using ADMAT direct way
xADobj = deriv(x,eye(length(x)));
ADfun = @(x)nlconstr(x,model,pvec);
xADres = ADfun(xADobj);
jac = sparse(getydot(xADres)); 

% method 2 - jacobian using ADMAT sparse jacobians
% Extra.fh = @(x)nlconstr(x,model,pvec);
% jpatt = getjpi('ADMATjaceg',3,6,Extra,'f');
% [~,jac] = evalj('ADMATjaceg',x,Extra,3,jpatt);


function jstr = nlconstrjacsparsity(model,pvec)
Extra.fh = @(x)nlconstr(x,model,pvec);
[~,jstr] = getjpi('ADMATjaceg',3,6,Extra,'f');
% Jpatt = modelSparsity(model);
% Jmc = Jpatt(1:3,1:3);
% Jpvec = sparse([1 0 1;0 1 1;0 0 0]);
% jspr = [Jmc Jpvec];



