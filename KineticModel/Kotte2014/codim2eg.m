% codim 2 bifurcation with CL_MATCONT 
% Example: catalytic oscillator
pvec = [2.5,1.92373,10,.0675,1,.1,.4];

% integrate to get initial equilibrium
odefun = @(t,x)oscillatorODE(t,x,[],pvec');
opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
tspan = 0:.1:500;
[tout,yout] = ode45(odefun,tspan,[0.00146723;0.826167;0.123119],opts);
figure
plot(tout,yout);
% [tout2,yout2] = ode45(odefun,tspan,[0.0029538;0.76211;0.16781],opts);
% figure
% plot(tout2,yout2);
% close all
ap = 2;
[data,y,p] = execMATCONT(@oscillator,[],yout(end,:)',pvec',ap,[],[],300);

% continue from limit points obtained on equilibrium branch
% collect only limit points
id = cat(1,data.s1.index);
label = cellstr(cat(1,data.s1.label));
label = cellfun(@(x)strcmpi(x,'LP'),label);
LPid = id(label);

LPx = data.x1(1:3,LPid(1));
LPp = p(LPid(1));
pvec(ap) = LPp;

global sys
sys.gui.pausespecial=1;  %Pause at special points 
sys.gui.pausenever=0;    %Pause never 
sys.gui.pauseeachpoint=0; %Pause at each point

LPap = [2 7];
[x0,v0]=init_LP_LP(@oscillator,LPx,pvec',LPap);
opt=contset;
opt=contset(opt,'MaxNumPoints',300);
opt=contset(opt,'MinStepSize',0.00001);
opt=contset(opt,'MaxStepSize',0.01);
opt=contset(opt,'Singularities',1);
opt = contset(opt,'Eigenvalues',1);
[x,v,s,h,f]=cont(@limitpoint,x0,v0,opt);
opt=contset(opt,'Backward',1);
[x2,v2,s2,h2,f2]=cont(@limitpoint,x0,v0,opt);

hfig = bifurcationPlot(x,s,f,[4,1]);

% calculate jacobian for points on continuation curve
[jac,lambda,w] = getKotteJacobian(@oscillatorNLAE,yout(end,:)',p',[])
