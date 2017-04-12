% inverse bifurcation from Lu, et al., 2006
% initial parameter sets for G1S Transition Model
k1 = 1;k2 = 1.6;k3 = .05;
k16 = .4;k61 = .3;k34 = .04;k43 = .01;k67 = .7;k76 = .1;
k23 = .3;k25 = .7;k28 = .06;
k89 = .07;k98 = .01;
a = .04;
j11 = .5;j12 = 5;j13 = .002;j15=.001;j18 = .6;
j61 = 5;j62 = 8;j63 = 2;j65 = 6;j68 = 7;
Km1 = .5;Km2 = 4;Km4 = .3;Km9 = .005;kp = .05;
phipRB = .005;phiE2F1 = .1;phiCycDi = .023;phiCycDa = .03;
phiAP1 = .01;phipRBp = .06;phipRBpp = .04;phiCycEi = .06;phiCycEa = .05;
Fm = 0;

pvec = [k1,k2,k3,k16,k61,k34,k43,k67,k76,k23,...
        k25,k28,k89,k98,a,j11,j12,j13,j15,j18,...
        j61,j62,j63,j65,j68,...
        Km1,Km2,Km4,Km9,kp,...
        phipRB,phiE2F1,phiCycDi,phiCycDa,phiAP1,...
        phipRBp,phipRBpp,phiCycEi,phiCycEa,Fm];

% parameter bounds    
% psl = [1e-4,0];
% psu = [1e-1,2];

% solve ode system to get initial equilibrium points
% ival = [0.001;.014;0;0;0;0;0;0;.001];
ival = [1;0;0;0;0;0];
% aefun = @(x)g1sNLAE(x,[],pvec);
% options = optimoptions(@fsolve,'Display','iter',...
%                         'TolFun',1e-15,'TolX',1e-15,...
%                         'MaxFunEvals',10000,'MaxIter',10000);
% [val,fval,eflag] = fsolve(aefun,ival,options);

nvar = 6;
tspan = 0:0.1:2000;
% Fmval = linspace(.001,0.5,100);
% allyout = zeros(nvar,length(Fmval));
ap = 40;

odefun = @(t,x)g1sODE(t,x,[],pvec');
opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
[tout,yout] = ode15s(odefun,tspan,ival,opts);
figure
plot(tout,yout);
legend(gca,'pRB','E2F1','CycDi','CycDa','AP1','pRBp'); 

% for iFm = 1:length(Fmval)
%     pvec(ap) = Fmval(iFm);
%     odefun = @(t,x)g1sODE(t,x,[],pvec');
%     opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
%     [tout,yout] = ode15s(odefun,tspan,ival,opts);
%     figure
%     plot(tout,yout);
%     legend(gca,'pRB','E2F1','CycDi','CycDa','AP1','pRBp'); 
%     allyout(:,iFm) = yout(end,:)';
%     ival = yout(end,:)';
%     drawnow
%     close all
% end
    
 % ,'pRBpp','CycEi','CycEa');
% [eigvals,eigw] = stabilityInfo(@g1sNLAE,yout(end,:),[],pvec');

% run MATCONT


% global sys
% sys.gui.pausespecial=0;  %Pause at special points 
% sys.gui.pausenever=1;    %Pause never 
% sys.gui.pauseeachpoint=0; %Pause at each point
% 
% % continuation from initial equilibrium - initialization
% % ap = 12; % index for parameter to be continued on     
% [x0,v0] = init_EP_EP(@g1s,yout(end,:)',pvec',ap);
% opt = contset;
% % opt = contset(opt,'VarTolerance',1e-4);
% % opt = contset(opt,'VarTolerance',1e-4);
% % opt = contset(opt,'FunTolerance',1e-4);
% opt = contset(opt,'MaxNumPoints',2000);
% opt = contset(opt,'MaxStepsize',0.01);
% opt = contset(opt,'MinStepSize',0.001);
% opt = contset(opt,'Singularities',1);
% opt = contset(opt,'Eigenvalues',1);
% opt = contset(opt,'Backward',0);
% 
% [x1,v1,s1,h1,f1] = cont(@equilibrium,x0,v0,opt); 




[data,y,p] = execMATCONT(@g1s,[],yout(end,:)',pvec',ap,[],[],3000);
figure
cpl(data.x1,data.v1,data.s1,[7 2]);
hold on
BPid = data.s1(2).index;

global sys
sys.gui.pausespecial=0;  %Pause at special points 
sys.gui.pausenever=1;    %Pause never 
sys.gui.pauseeachpoint=0; %Pause at each point

pvec(ap) = data.x1(end,BPid);
xBP = data.x1(1:end-1,BPid);
opt = contset;
opt = contset(opt,'VarTolerance',1e-4);
opt = contset(opt,'VarTolerance',1e-4);
opt = contset(opt,'FunTolerance',1e-4);
opt = contset(opt,'MaxNumPoints',1000);
opt = contset(opt,'MaxStepsize',0.01);
opt = contset(opt,'MinStepSize',0.00001);
opt = contset(opt,'Singularities',1);
opt = contset(opt,'Eigenvalues',1);
opt=contset(opt,'Backward',1);

[x0,v0] = init_BP_BP(@g1s,xBP,pvec',[1 11 30],ap);
[x3,v3,s3,h3,f3]=cont(@branchpoint,x0,[],opt);


% perform continuation from Hopf bifurcation point on equilibrium cont curve
apH = [2 1];

% extract all limit points from data
id = cat(1,data.s1.index);
label = cellstr(cat(1,data.s1.label));
Hid = id(cellfun(@(x)strcmpi(x,'H'),label));

% global sys
% sys.gui.pausespecial=1;  %Pause at special points 
% sys.gui.pausenever=0;    %Pause never 
% sys.gui.pauseeachpoint=0; %Pause at each point

% pvec(ap) = data.x1(end,Hid(1));
% [x0,v0]=init_H_H(@repressilator,data.x1(1:6,Hid(1)),pvec',apH);
% % set continuation options
% opt=contset;
% % opt = contset(opt,'VarTolerance',1e-3);
% % opt = contset(opt,'VarTolerance',1e-3);
% % opt = contset(opt,'FunTolerance',1e-3);
% opt=contset(opt,'MaxNumPoints',100000);
% opt=contset(opt,'MinStepSize',0.0001);
% opt=contset(opt,'MaxStepSize',0.1);
% opt=contset(opt,'Singularities',1);
% opt = contset(opt,'Eigenvalues',1);
% % opt = contset(opt,'Adapt',1);
% % opt = contset(opt,'Multipliers',1);
% [x,v,s,h,f]=cont(@hopf,x0,v0,opt);

% limit cycle conitnuation (1 free variables in apH)
% apH = 2;
% pvec(ap) = data.x1(end,Hid);
% [x0,v0]=init_H_LC(@repressilator,data.x1(1:6,Hid),pvec',apH,1e-8,30,4);% 
% % set continuation options
% opt=contset;
% opt = contset(opt,'VarTolerance',1e-3);
% opt = contset(opt,'VarTolerance',1e-3);
% opt = contset(opt,'FunTolerance',1e-3);
% opt=contset(opt,'MaxNumPoints',300);
% opt=contset(opt,'MinStepSize',0.001);
% opt=contset(opt,'MaxStepSize',10);
% opt=contset(opt,'Singularities',1);
% opt = contset(opt,'Eigenvalues',1);
% opt = contset(opt,'Adapt',1);
% opt = contset(opt,'Multipliers',1);
% [x,v,s,h,f]=cont(@limitcycle,x0,v0,opt);


%  Hdata =...
%  execHcont(@oscillator,yout(end,:)',pvec',apH,ap,data);
%  bifurcationPlot(Hdata.x1,Hdata.s1,Hdata.f1,[4,1],[],1,fig);

