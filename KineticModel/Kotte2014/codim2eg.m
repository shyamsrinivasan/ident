% codim 2 bifurcation with CL_MATCONT 
% Example: catalytic oscillator
p = [2.5,1.92373,10,.0675,1,.1,.4];

% integrate to get initial equilibrium
odefun = @(t,x)oscillatorODE(t,x,p');
opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
tspan = 0:.1:500;
[tout,yout] = ode45(odefun,tspan,[0.00146723;0.826167;0.123119],opts);
figure
plot(tout,yout);
% [tout2,yout2] = ode45(odefun,tspan,[0.0029538;0.76211;0.16781],opts);
% figure
% plot(tout2,yout2);
% close all

global sys
sys.gui.pausespecial=1;  %Pause at special points 
sys.gui.pausenever=0;    %Pause never 
sys.gui.pauseeachpoint=0; %Pause at each point

ap = 2;
[x0,v0]=init_EP_EP(@oscillator,yout(end,:)',p',ap);
opt=contset;
opt=contset(opt,'MaxNumPoints',300);
opt=contset(opt,'MinStepSize',0.00001);
opt=contset(opt,'MaxStepSize',0.01);
opt=contset(opt,'Singularities',1);
% opt = contset(opt,'Eigenvalues',1);
[x,v,s,h,f]=cont(@equilibrium,x0,v0,opt);
