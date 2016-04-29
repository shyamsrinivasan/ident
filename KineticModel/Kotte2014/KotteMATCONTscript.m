% Kotte MATCONT
% run Kotte2014_script before running this script

% continuation and dynamical systems analysis using MATCONT

% remove metabolites held constant from consideration in the model
% integration phase
[model,pvec,newmc,cnstmet] =...
remove_eMets(FBAmodel,parameter,mc,[FBAmodel.Vind FBAmodel.Vex],...
{'enz1[c]','enz1[e]','enz[e]','ac[e]','bm[c]','bm[e]','pep[e]'});

% only initialize for varmets   
nvar = length(model.mets)-length(find(cnstmet));
M = newmc(1:nvar);
PM = newmc(nvar+1:end);
model.PM = PM;

% parameters
clear pvec
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
    
% Kotte_givenscript
allhandles = feval(@Kotte2014glycolysis);
rhsfunc = allhandles{2};
givenModel = @(t,x)rhsfunc(t,x,model,kEcat,KEacetate,...
        KFbpFBP,vFbpmax,Lfbp,KFbpPEP,...
        vEXmax,KEXPEP,...
        vemax,KeFBP,ne,acetate,d);
fluxg = Kotte_givenFlux([M;model.PM],pvec,model);
dMdtg = givenModel(0,M);

% given model SS
opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
[tout,yout] = ode45(givenModel,0:0.1:200,M,opts);
fout = zeros(length(tout),4);
for it = 1:length(tout)
    fout(it,:) = Kotte_givenFlux([yout(it,:)';model.PM],pvec,model);
end
xeq = yout(end,:)';

% given model fsolve
gfun = @(x)Kotte_givenNLAE(x,model,pvec);
dMg = gfun(M);
options = optimoptions('fsolve','Display','iter','TolFun',1e-10,'TolX',1e-10);
[x1,fval,exitflag,output,jacobian] = fsolve(gfun,M,options);
fgout = Kotte_givenFlux([x1;model.PM],pvec,model);

% continuation and dynamical systems analysis using MATCONT
global sys
sys.gui.pausespecial=0;  %Pause at special points 
sys.gui.pausenever=1;    %Pause never 
sys.gui.pauseeachpoint=0; %Pause at each point

% continuation from initial equilibrium - initialization
ap = 12; % index for parameter to be continued on     
[x0,v0] = init_EP_EP(@KotteMATCONT,xeq,pvec,ap);

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

% flux calculation
flux1 = zeros(4,size(x1,2));
y = x1(1:length(xeq),:);
p = x1(length(xeq)+1:end,:);
for it = 1:size(x1,2)
    pvec(ap) = p(it);
    flux1(:,it) = KotteMATCONTflux(y(:,it),pvec);
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

% figure
% bifurcationPlot(flux1,flux1,s1,f1,2,1);