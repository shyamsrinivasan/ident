M = zeros(3,1);
M(1)  = 0.001;   % PEP
M(2)  = 10;   % FBP
M(3)  = 1;      % E

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
KeFBP = 0.45;       % or 0.45
ne = 2;             % or 2
acetate = 0.1;      % a.u acetate
d = 0.25;           % or 0.25 or 0.35
pvec = [kEcat,KEacetate,...
        KFbpFBP,vFbpmax,Lfbp,KFbpPEP,...
        vEXmax,KEXPEP,...
        vemax,KeFBP,ne,acetate,d];

allhandles = feval(@Kotte2014glycolysis);
rhsfunc = allhandles{2};
func = @(t,x)rhsfunc(t,x,kEcat,KEacetate,...
        KFbpFBP,vFbpmax,Lfbp,KFbpPEP,...
        vEXmax,KEXPEP,...
        vemax,KeFBP,ne,acetate,d);
dMdt = func(0,M);    

opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
[tout,yout] = ode45(func,0:0.1:200,M,opts);
fout = zeros(length(tout),4);
for it = 1:length(tout)
    fout(it,:) = Kotte_glycolysisflux(yout(it,:),pvec);
end
plotKotteVariables(tout,yout,1);
plotKotteVariables(tout,fout,2);
xeq = yout(end,:)';
close all

% perturbation to equilibrium solution
% positive perturbations
% t = linspace(0,60,1000);
% allYout = zeros(length(t),3,1000);
% Nxeq = eqsolperturbation(xeq,n,'pos');
% Nxeq = repmat(xeq,1,1000)+...
%        repmat(xeq,1,1000).*random(makedist('Uniform'),length(xeq),1000);
% for i = 1:1000    
%     [tout,yout] = ode45(func,t,Nxeq(:,i));
%     allYout(:,:,i) = yout;    
% end

% negative perturbations
% t = linspace(0,60,1000);
% allYout = zeros(length(t),3,1000);
% Nxeq = eqsolperturbation(xeq,n,'neg');
% Pxeq = zeros(3,1000);
% Nxeq = repmat(xeq,1,1000)-...
%        repmat(xeq,1,1000).*random(makedist('Uniform'),length(xeq),1000);
% for i = 1:1000    
%     Pxeq(:,i) = max(0.00001,Nxeq(:,i));
%     [tout,yout] = ode45(func,t,Pxeq(:,i));
%     allYout(:,:,i) = yout;    
% end

%sign randomized perturbations
n = 100000;
% t = linspace(0,60,1000);
% allYout = zeros(length(t),3,n);
% Nxeq = eqsolperturbation(xeq,n,'rnd');
% Pxeq = zeros(3,n);
% Nxeq = repmat(xeq,1,n)+...
%        randi([-1 1],3,n).*repmat(xeq,1,n).*(2*random(makedist('Uniform'),length(xeq),n));
% for i = 1:n    
%     Pxeq(:,i) = max(0.00001,Nxeq(:,i));
%     [tout,yout] = ode45(func,t,Pxeq(:,i));
%     allYout(:,:,i) = yout;    
% end
% plotKotteVariables(tout,allYout,1);

% NL solve for rhs of Kotte ode
% fun = @(x)Kotte_glycolysis_NLAE(x,pvec);
% options = optimoptions('fsolve','Display','iter','TolFun',1e-10,'TolX',1e-10);
% [x1,fval,exitflag,output,jacobian] = fsolve(fun,xeq,options);
% flux = Kotte_glycolysisflux(x1,pvec);
% display(x1);
% display(flux);
% 
% [tout,yout] = ode45(func,0:0.1:1000,x1,opts);
% fout = zeros(length(tout),4);
% for it = 1:length(tout)
%     fout(it,:) = Kotte_glycolysisflux(yout(it,:),pvec);
% end
% plotKotteVariables(tout,yout,1);
% plotKotteVariables(tout,fout,2);
% close all
% 
% 
% change in initial guess and acetate concentration for fsolve
% lb = [.001;0.001;0.001];
% ub = [1;10;10];
% Pxeq = zeros(3,n);
% acetate = 0.001 + (0.2-0.001).*random(makedist('Uniform'),1,n);
% acetate = sort(acetate);
% uniqueAcetate = acetate(1);
% for jac = 1:length(acetate)
%     idAcetate = abs(uniqueAcetate-repmat(acetate(jac),size(uniqueAcetate,1),1))<1e-4;
%     if ~any(idAcetate)
%         uniqueAcetate = [uniqueAcetate;acetate(jac)];
%     else
%         acetate(jac) = uniqueAcetate(logical(idAcetate));
%     end
% end
% acetate = acetate';
% 
% xNLEout = zeros(n,3);
% fNLEout = zeros(n,4);
% allFlag = zeros(n,1);
% Nxeq = repmat(lb,1,n) + (repmat(ub,1,n)-repmat(lb,1,n)).*...
%        random(makedist('Uniform'),length(xeq),n);
% Nxeq = ivalueperturbation(lb,ub,n);
% options = optimoptions(options,'Display','off');   
% for i = 1:n
%     fprintf('Iteration : %d\n\n',i);
%     Pxeq(:,i) = max(0.00001,Nxeq(:,i));
%     pvec(12) = acetate(i);    
%     [xNLE,fval,exitflag,output,jacobian] = fsolve(fun,Pxeq(:,i),options);
%     flux = Kotte_glycolysisflux(xNLE,pvec);
%     xNLEout(i,:) = xNLE;
%     fNLEout(i,:) = flux;
%     allFlag(i) = exitflag;
% end   
% xNLEout(allFlag<=0,:) = [];
% fNLEout(allFlag<=0,:) = [];
% Pxeq(:,allFlag<=0) = [];
% acetate(allFlag<=0) = [];
% allFlag(allFlag<=0) = [];
% 
% fNLEout(min(xNLEout,[],2)<0,:) = [];
% Pxeq(:,min(xNLEout,[],2)<0) = [];
% allFlag(min(xNLEout,[],2)<0) = [];
% acetate(min(xNLEout,[],2)<0) = [];
% xNLEout(min(xNLEout,[],2)<0,:) = [];
% 
% newVector = zeros(size(xNLEout,1),size(xNLEout,2)+size(Pxeq,1)+1+size(fNLEout,1));
% uniqueSS = xNLEout(1,:);
% for j = 1:size(xNLEout,1)  
%     ind = prod(abs(uniqueSS-repmat(xNLEout(j,:),size(uniqueSS,1),1))<1e-4,2);
%     if ~any(ind)
%         % add to uniqueSS as a new ss and do not replace
%         uniqueSS = [uniqueSS;xNLEout(j,:)];
%     else
%         % subsitute with match in uniqueSS        
%         xNLEout(j,:) =...
%         uniqueSS(logical(ind),:);
%         fNLEout(j,:) = Kotte_glycolysisflux(uniqueSS(logical(ind),:),pvec);        
%     end
%     newVector(j,:) = [xNLEout(j,:) Pxeq(:,j)' acetate(j) fNLEout(j,:)];
% end
        
% plotKotteVariables(Pxeq,allNLy,3);  


% change in acetate concentration
% acetate = linspace(0.001,3,n);
% xNLEout = zeros(n,3);
% fNLEout = zeros(n,4);
% allFlag = zeros(n,1);
% for i = 1:n
%     pvec(12) = acetate(i);    
%     [xNLE,fval,exitflag,output,jacobian] = fsolve(fun,xeq,options);
%     flux = Kotte_glycolysisflux(xNLE,pvec);
%     xNLEout(i,:) = xNLE;
%     fNLEout(i,:) = flux;
%     allFlag(i) = exitflag;
% end
% figure
% plot(acetate,fNLEout(:,1));

% fsolve from random initial conditions
% lb = [.001;0.001;0.001];
% ub = [1;10;10];
% Pxeq = zeros(3,n);
% acetate = linspace(2,0.001,n);
% xNLEout = zeros(n,3);
% fNLEout = zeros(n,4);
% allFlag = zeros(n,1);
% Nxeq = repmat(lb,1,n) + (repmat(ub,1,n)-repmat(lb,1,n)).*...
%        random(makedist('Uniform'),length(xeq),n);
% for i = 1:n
%     fprintf('Iteration : %d\n\n',i);
%     Pxeq(:,i) = max(0.00001,Nxeq(:,i));
%     pvec(12) = acetate(i);    
%     [xNLE,fval,exitflag,output,jacobian] = fsolve(fun,Pxeq(:,i),options);
%     flux = Kotte_glycolysisflux(xNLE,pvec);
%     xNLEout(i,:) = xNLE;
%     fNLEout(i,:) = flux;
%     allFlag(i) = exitflag;
% end   
% xNLEout(allFlag<=0,:) = [];
% fNLEout(allFlag<=0,:) = [];
% Pxeq(:,allFlag<=0) = [];
% allFlag(allFlag<=0) = [];
% 
% fNLEout(min(xNLEout,[],2)<0,:) = [];
% Pxeq(:,min(xNLEout,[],2)<0) = [];
% allFlag(min(xNLEout,[],2)<0) = [];
% xNLEout(min(xNLEout,[],2)<0,:) = [];

% continuation and dynamical systems analysis using MATCONT
global sys
sys.gui.pausespecial=0;  %Pause at special points 
sys.gui.pausenever=1;    %Pause never 
sys.gui.pauseeachpoint=0; %Pause at each point

% continuation from initial equilibrium - initialization
ap = 12; % index for parameter to be continued on     
[x0,v0] = init_EP_EP(@Kotte2014glycolysis,xeq,pvec,ap);

% MATCONT options
opt = contset;
opt = contset(opt,'VarTolerance',1e-3);
opt = contset(opt,'VarTolerance',1e-3);
opt = contset(opt,'FunTolerance',1e-3);
opt = contset(opt,'MaxNumPoints',500);
opt = contset(opt,'MaxStepsize',.01);
opt = contset(opt,'Singularities',1);
opt = contset(opt,'Eigenvalues',1);

% Equilibrium Continuation
[x1,v1,s1,h1,f1] = cont(@equilibrium,x0,v0,opt); 
figure
% flux calculation
flux1 = zeros(4,size(x1,2));
y = x1(1:length(xeq),:);
p = x1(length(xeq)+1:end,:);
for it = 1:size(x1,2)
    pvec(ap) = p(it);
    flux1(:,it) = Kotte_glycolysisflux(y(:,it),pvec);
end

figure
subplot(221)
bifurcationPlot(y,p,s1,f1,1,1)
xlabel('Acetate');
ylabel('E');
subplot(222)
bifurcationPlot(y,p,s1,f1,2,1)
xlabel('Acetate');
ylabel('PEP');
subplot(223)
bifurcationPlot(y,p,s1,f1,3,1)
xlabel('Acetate');
ylabel('FBP');
subplot(224)
bifurcationPlot(flux1,p,s1,f1,1,1);
xlabel('Acetate');
ylabel('flux J');

figure
subplot(221)
bifurcationPlot(flux1,p,s1,f1,4,1);
xlabel('Acetate');
ylabel('flux E');
subplot(222)
bifurcationPlot(y,flux1,s1,f1,1,1);
xlabel('flux J');
ylabel('E');
subplot(223)
bifurcationPlot(y,flux1,s1,f1,2,1);
xlabel('flux J');
ylabel('PEP');
subplot(224)
bifurcationPlot(y,flux1,s1,f1,3,1);
xlabel('flux J');
ylabel('FBA');


% figure
% cpl(x1,v1,s1,[4,1]);


% extract critical points - if critical points exist
ap1 = 9; % index for second parameter to be continued on for a pitch fork LP
if size(s1,1)>2
    cindx = cat(1,s1.index);
    cindx = cindx(2:end-1);
    for inc = 1:length(cindx)
        opt=contset(opt,'Backward',[]);
        cid = cindx(inc);
        xCP = x1(1:length(xeq),cid);
        pvec(ap) = x1(length(xeq)+1:end,cid);
        newap = [ap ap1];
        [x0,v0] = init_LP_LP(@Kotte2014glycolysis,xCP,pvec',newap);
        [x2,v2,s2,h2,f2]=cont(@limitpoint,x0,v0,opt);
        hold on
        cpl(x2,v2,s2,[4,1]);
        opt=contset(opt,'Backward',1);
        [x3,v3,s3,h3,f3]=cont(@limitpoint,x0,v0,opt);
        cpl(x3,v3,s3,[4,1]);
    end
end
    
% 
% %extreact branch point
% xBP = x1(1:length(xeq),s1(2).index);
% pvec(ap) = x1(4,s1(2).index);
% 
% [x0,v0] = init_BP_EP(@Kotte2014glycolysis,xBP,pvec,s1(2),0.01);
% [x2,v2,s2,h2,f2]=cont(@equilibrium,x0,v0,opt); %Switch branches and continue.
% figure()
% cpl(x2,v2,s2,[2 1]);


