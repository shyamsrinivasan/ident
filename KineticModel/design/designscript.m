% inverse bifurcation from Lu, et al., 2006
% initial parameter sets
alpha = 100; % 10^2.5;
beta = 1;
delta = 1e-3;
h = 1.5;
pvec = [alpha,beta,delta,h];
pi = [alpha,beta];
ps = [delta,h];

% parameter bounds    
psl = [1e-4,0];
psu = [1e-1,2];

% solve ode system to get initial equilibrium points
ival = [10;10;10;10;10;10];
tspan = 0:0.1:100;
odefun = @(t,x)repressilatorODE(t,x,[],pvec');
opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
[tout,yout] = ode15s(odefun,tspan,ival,opts);
figure
plot(tout,yout);

% run MATCONT
ap = 2;
[data,y,p] = execMATCONT(@repressilator,[],yout(end,:)',pvec',ap,[],[],300);

% perform continuation from Hopf bifurcation point on equilibrium cont curve


