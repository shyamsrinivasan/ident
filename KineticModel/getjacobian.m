function [J,lambda,w] = getjacobian(x,pvec,model)
% get analytical jacobian based on modular rate law
nrxn = model.nt_rxn;
nmet = length(x);
DVX = zeros(nmet,nrxn);

Vind = model.Vind;
Vex = model.Vex;

DVX(:,Vind) = CKjacobian(model,pvec,x,Vind);
DVX(:,Vex) = TKjacobian(model,pvec,x,Vex);

S = model.S;
J = sparse(zeros(nmet,nmet));

% complete jacobian - build in column loops rowwise
% J = nmet x nmet
for imet = 1:nmet
    J(:,imet) = S*DVX(imet,:)';
end

% since extracellular concentrations are sometime zero the jacobian for
% internal metabolites alone is also given separately
nint = model.nint_metab;
Jint = J(1:nint,1:nint);

% use ADMAT to calculate jacobians
admatfun = @(x)toyNLAE(x,model,pvec);
xADMATobj = deriv(x,eye(size(x,1)));
xADMATres = admatfun(xADMATobj);
F = getval(xADMATres);
J = getydot(xADMATres); 

% eigen value and eigen vector
[w,lambda] = eig(J);
lambda = diag(lambda);
