% optimization framework to estimate mss producing/non-producing enzyme
% parameters

cplexlp
Opt = opti('fun',fun,'bounds',lb,ub,'x0',x0);

Opt = opti('f',f,'ineq',A,b,'bounds',lb,ub);
% f = c;
% A = [a11 a12;a21 a22];
% b = [0;0];
% lb = [0;0];
% ub = [0;0];
[x,fval,exitflag,info] = solve(Opt);

% fun - nonlinear objective function
Opt = opti('f',f,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'nljacstr',nljacstr);

[x,fval,exitflag,info] = solve(Opt,x0);

% use ADMAT to calculate jacobians
fun = ADfun('',length(M));
y = deriv(M,diag(ones(3,1)));

