function jac = HoltJacobian(y)

nvar = size(y,1);
ADMATobj = deriv(y,eye(nvar));
ADMATres = Holt(ADMATobj);
jac = getydot(ADMATres);