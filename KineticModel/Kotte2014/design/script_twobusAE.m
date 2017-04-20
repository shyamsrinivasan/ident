% run twobusAE to find equilibrium points and get fold bifurcation manifold
x0 = [-0.1;-1];
p = [.1;.1];

opts = odeset('AbsTol',1e-8,'RelTol',1e-10);
odefun = @(t,x)twobusODE(t,x,p);
[t,y] = ode45(odefun,0:.1:1,x0);
figure
plot(t,y);


