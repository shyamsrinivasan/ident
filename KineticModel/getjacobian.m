function [J,lambda,w] = getjacobian(x,pvec,model)
% use ADMAT to calculate jacobians

admatfun = @(x)toyNLAE(x,model,pvec);
xADMATobj = deriv(x,eye(size(x,1)));
xADMATres = admatfun(xADMATobj);
F = getval(xADMATres);
J = getydot(xADMATres); 

% eigen value and eigen vector
[w,lambda] = eig(J);
lambda = diag(lambda);
