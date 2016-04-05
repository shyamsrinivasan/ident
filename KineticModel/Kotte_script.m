M = zeros(3,1);
M(1)  = 1;      % E
M(2)  = 0.01;   % PEP
M(3)  = 0.03;   % FBP

%parameters
kEcat = 1;
KEacetate = 0.1;    % or 0.02
KFbpFBP = 0.1;
vFbpmax = 1;
Lfbp = 4e6;
KFbpPEP = 0.1;
vEXmax = 1;
KEXPEP = 0.3;
vemax = 1.1;        % for bifurcation analysis: 0.7:0.1:1.3
KeFBP = 0.1;        % or 0.45
ne = 1;             % or 2
acetate = 0.1;        % a.u acetate
d = 0.18;           % or 0.25 or 0.35
pvec = [kEcat,KEacetate,...
        KFbpFBP,vFbpmax,Lfbp,KFbpPEP,...
        vEXmax,KEXPEP,...
        vemax,KeFBP,ne,acetate,d];

dMdt = Kotte_glycolysis(0,M,pvec);


% Xss = J*KEXPEP/(vEXmax-J);
func = @(t,x)Kotte_glycolysis(t,x,pvec);
[tout,yout] = ode45(func,0:0.1:60,M);
subplot(221);
plot(tout,yout(:,2));
ylabel('PEP');
hold on
subplot(222);
plot(tout,yout(:,3));
ylabel('FBP');
hold on
subplot(223);
plot(tout,yout(:,1));
ylabel('super Enzyme E');
xlabel('time');
hold on

xeq = yout(end,:)';
% perturbation to equilibrium solution
% positive perturbations
% t = linspace(0,60,1000);
% allYout = zeros(length(t),3,1000);
% Nxeq = repmat(xeq,1,1000)+...
%        repmat(xeq,1,1000).*random(makedist('Uniform'),length(xeq),1000);
% for i = 1:1000    
%     [tout,yout] = ode45(func,t,Nxeq(:,i));
%     allYout(:,:,i) = yout;    
% end

% negative perturbations
% t = linspace(0,60,1000);
% allYout = zeros(length(t),3,1000);
% Pxeq = zeros(3,1000);
% Nxeq = repmat(xeq,1,1000)-...
%        repmat(xeq,1,1000).*random(makedist('Uniform'),length(xeq),1000);
% for i = 1:1000    
%     Pxeq(:,i) = max(0.00001,Nxeq(:,i));
%     [tout,yout] = ode45(func,t,Pxeq(:,i));
%     allYout(:,:,i) = yout;    
% end

%sign randomized perturbations
n = 5000;
% t = linspace(0,60,1000);
% allYout = zeros(length(t),3,n);
% Pxeq = zeros(3,n);
% Nxeq = repmat(xeq,1,n)+...
%        randi([-1 1],3,n).*repmat(xeq,1,n).*(2*random(makedist('Uniform'),length(xeq),n));
% for i = 1:n    
%     Pxeq(:,i) = max(0.00001,Nxeq(:,i));
%     [tout,yout] = ode45(func,t,Pxeq(:,i));
%     allYout(:,:,i) = yout;    
% end
% 
% 
% for i = 1:1000
%     subplot(221);
%     plot(tout,allYout(:,2,i),'r');
%     ylabel('PEP');
%     hold on
%     subplot(222);
%     plot(tout,allYout(:,3,i),'r');
%     ylabel('FBP');
%     hold on
%     subplot(223);
%     plot(tout,allYout(:,1,i),'r');
%     ylabel('super Enzyme E');
%     xlabel('time');
%     hold on
% end
% flux calculation
% flux = zeros(length(tout),4);
% for it = 1:length(tout)
%     flux(it,:) = Kotte_glycolysisflux(yout(it,:),pvec);
% end

% subplot(221);
% plot(tout,yout(:,2),'r');
% ylabel('PEP');
% hold on
% subplot(222);
% plot(tout,yout(:,3),'r');
% ylabel('FBP');
% hold on
% subplot(223);
% plot(tout,yout(:,1),'r');
% ylabel('super Enzyme E');
% xlabel('time');
% hold on

%phase planes
% figure
% subplot(131);
% plot(yout(:,1),yout(:,2));
% xlabel('super Enzyme E');
% ylabel('PEP');
% subplot(132);
% plot(yout(:,1),yout(:,3));
% xlabel('super Enzyme E');
% ylabel('FBP');
% subplot(133);
% plot(yout(:,2),yout(:,3));
% xlabel('PEP');
% ylabel('FBP');

% figure
% subplot(221);
% plot(tout,flux(:,2));
% ylabel('Fbp');
% subplot(222);
% plot(tout,flux(:,3));
% ylabel('FBP');
% subplot(223);
% plot(tout,flux(:,4));
% ylabel('E(FBP)');
% xlabel('time');
% subplot(224);
% plot(tout,flux(:,1));
% ylabel('super Enzyme E,J');
% xlabel('time');
% close all

%NL solve for rhs of Kotte ode
fun = @(x)Kotte_glycolysis_NLAE(x,pvec);
options = optimoptions('fsolve','Display','iter');
[x1,fval,exitflag,output,jacobian] = fsolve(fun,M,options);

%change in acetate concentration
acetate = linspace(0.001,100,10000);
xNLEout = zeros(10000,3);
fNLEout = zeros(10000,4);
for i = 1:10000
    pvec(12) = acetate(i);    
    [xNLE,fval,exitflag,output,jacobian] = fsolve(fun,M,options);
    flux = Kotte_glycolysisflux(M,pvec);
    xNLEout(i,:) = xNLE;
    fNLEout(i,:) = flux;
end
% figure
% plot(acetate,xNLEout(:,1));

%change in initial guess for fsolve
% allNLy = zeros(3,n);
% Pxeq = zeros(3,n);
% Nxeq = repmat(xeq,1,n)+...
%        randi([-1 1],3,n).*repmat(xeq,1,n).*random(makedist('Uniform'),length(xeq),n);
% options = optimoptions('fsolve','Display','off');
% for i = 1:n
%     fprintf('Iteration : %d\n\n',i);
%     Pxeq(:,i) = max(0.00001,Nxeq(:,i));
%     [xNLE,fval,exitflag,output,jacobian] = fsolve(fun,Pxeq(:,i),options);
%     allNLy(:,i) = xNLE;    
% end


% figure
% subplot(221);
% plot(Pxeq(2,:),allNLy(2,:),'LineStyle','none','Marker','o','MarkerSize',5);
% ylabel('PEP');
% subplot(222);
% plot(Pxeq(3,:),allNLy(3,:),'LineStyle','none','Marker','o','MarkerSize',5);
% ylabel('FBP');
% subplot(223);
% plot(Pxeq(1,:),allNLy(1,:),'LineStyle','none','Marker','o','MarkerSize',5);
% ylabel('super Enzyme E');



    

%continuation and dynamical systems analysis using MATCONT
%initial equilibrium
% ap = 9; % index for parameter to be continued on     
% [x0,v0] = init_EP_EP(@Kotte2014glycolysis,xeq',pvec,ap);
% 
% opt = contset;
% opt = contset(opt,'VarTolerance',1e-3);
% opt=contset(opt,'VarTolerance',1e-3);
% opt=contset(opt,'FunTolerance',1e-3);
% opt=contset(opt,'MaxNumPoints',500);
% opt=contset(opt,'MaxStepsize',.01);
% opt=contset(opt,'Singularities',1);
% opt=contset(opt,'Eigenvalues',1);
% [x1,v1,s1,h1,f1]=cont(@equilibrium,x0,[],opt); %Equilibrium continuation
% figure
% cpl(x1,v1,s1,[2,3,4]);
% 
% %extreact branch point
% xBP = x1(1:length(xeq),s1(2).index);
% pvec(ap) = x1(4,s1(2).index);
% 
% [x0,v0] = init_BP_EP(@Kotte2014glycolysis,xBP,pvec,s1(2),0.01);
% [x2,v2,s2,h2,f2]=cont(@equilibrium,x0,v0,opt); %Switch branches and continue.
% figure()
% cpl(x2,v2,s2,[2 1]);


