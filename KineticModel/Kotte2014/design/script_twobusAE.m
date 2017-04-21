% run twobusAE to find equilibrium points and get fold bifurcation manifold
x0 = [-0.18;.1];
p = [.5;.3];

opts = odeset('AbsTol',1e-8,'RelTol',1e-10);
odefun = @(t,x)twobusODE(t,x,p);
[t,y] = ode45(odefun,0:.1:10,x0,opts);
figure
plot(t,y);

% get stability info
jac1 = twobus_jacobian(x0,p);
[v1,eigval1,w1] = eig(jac1);
eigval1 = diag(eigval1);

jac2 = twobus_jacobian(y(end,:),p);
[v2,eigval2,w2] = eig(jac2);
eigval2 = diag(eigval2);

% get steady state using fsolve - gives a different equilibrium solution
% aefun = @(x)twobusAE(x,p);
% clear opts
% opts = optimoptions(@fsolve,'Display','iter','TolFun',1e-10,'TolX',1e-8);
% [x,fval,flag] = fsolve(aefun,x0,opts);


