function [J,lambda,w] = getKotteJacobian(funhdle,x,pvec,model)

% use ADMAT to calculate jacobians
admatfun = @(x)funhdle(x,model,pvec);

xADMATobj = deriv(x,eye(3));
xADMATres = admatfun(xADMATobj);
% F = getval(xADMATres);
J = getydot(xADMATres); 

% eigen value and eigen vector
[w,lambda] = eig(J);
lambda = diag(lambda);

