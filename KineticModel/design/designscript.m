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
apH = [2 1];

% extract all limit points from data
id = cat(1,data.s1.index);
label = cellstr(cat(1,data.s1.label));
Hid = id(cellfun(@(x)strcmpi(x,'H'),label));

% set continuation options
opt=contset;
opt=contset(opt,'MaxNumPoints',contpts);
opt=contset(opt,'MinStepSize',0.00001);
opt=contset(opt,'MaxStepSize',0.1);
opt=contset(opt,'Singularities',1);
opt = contset(opt,'Eigenvalues',1);

global sys
sys.gui.pausespecial=1;  %Pause at special points 
sys.gui.pausenever=0;    %Pause never 
sys.gui.pauseeachpoint=0; %Pause at each point

% limit point conitnuation (2 free variables in apH)
[x0,v0]=init_H_H(funame,x0,pvec,apH);
[x,v,s,h,f]=cont(@hopf,x0,v0,opt);


%  Hdata =...
%  execHcont(@oscillator,yout(end,:)',pvec',apH,ap,data);
%  bifurcationPlot(Hdata.x1,Hdata.s1,Hdata.f1,[4,1],[],1,fig);

