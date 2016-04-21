addpath(genpath('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel'));
rxfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\Kotte2014\Kotte2014.txt';
cnfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\Kotte2014\Kotte2014C.txt';

% create model structure
[FBAmodel,parameter,variable,nrxn,nmetab] = modelgen(rxfname);

% obtain conentrations from file
[mc,FBAmodel,met] = readCNCfromFile(cnfname,FBAmodel);

% remove metabolites held constant from consideration in the model
% integration phase
[model,pvec,newmc,cnstmet] =...
remove_eMets(FBAmodel,parameter,mc,[FBAmodel.Vind FBAmodel.Vex],...
{'enz1[c]','ac[e]','fdp[e]'});

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

% given model fsolve
gfun = @(x)Kotte_givenNLAE(x,model,pvec);
dMg = gfun(M);
options = optimoptions('fsolve','Display','iter','TolFun',1e-10,'TolX',1e-10);
[x1,fval,exitflag,output,jacobian] = fsolve(gfun,M,options);
fgout = Kotte_givenFlux([x1;model.PM],pvec,model);

% Kotte_Cscript
allhandles = feval(@Kotte2014Ckinetics);
rhsfunc = allhandles{2};
CModel = @(t,x)rhsfunc(t,x,model,kEcat,KEacetate,...
        KFbpFBP,vFbpmax,Lfbp,KFbpPEP,...
        vEXmax,KEXPEP,...
        vemax,KeFBP,ne,acetate,d);
fluxC = Kotte_CFlux([x1;model.PM],pvec,model);
dMdtC = CModel(0,x1);

% C model fsolve
cfun = @(x)Kotte_CNLAE(x,model,pvec);
dMc = cfun(x1);
options = optimoptions('fsolve','Display','iter',...
                       'TolFun',1e-10,...
                       'TolX',1e-10,...
                       'MaxFunEvals',1000000,...
                       'MaxIter',50000);
[x2,fval,exitflag,output,jacobian] = fsolve(cfun,x1,options);
fcout = Kotte_CFlux([M;model.PM],pvec,model);    

% C model SS
opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
[tout,yout] = ode45(CModel,0:0.1:60,M,opts);
fout = zeros(length(tout),4);
for it = 1:length(tout)
    fout(it,:) = Kotte_givenFlux([yout(it,:)';model.PM],pvec,model);
end

% pvec.d = 0.25;
% allhandles = feval(@Kotte2014Ckinetics);
% rhsfunc = allhandles{2};
% func = @(t,x)KotteCkinetics(t,x,model,pvec);
% dMdt = func(0,M);
allmc = [M;model.PM];

fun = @(x)KotteCkinetics(x,pvec,model);
dMdt = fun(M);


[x1,fval,exitflag,output,jacobian] = fsolve(fun,M,options);

% opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
% [tout,yout] = ode45(func,0:0.1:60,M,opts);
% opts = CVodeSetOptions('RelTol',1e-12,'AbsTol',1e-10);
% CVodeInit(func,'BDF','Newton',0.0,M,opts);
% [status,tout,yout] = CVode(0.1:0.1:15,'Normal');  
% plot(tout,yout);
% % calculate fluxes
% fout = zeros(length(tout),4);
% for it = 1:length(tout)
%     allmc = [yout(:,it);model.PM];
%     fout(it,[1;3;4]) = CKinetics(model,pvec,allmc,[1 3 4]);
%     fout(it,strcmpi(model.rxns,'ACpts')) = allmc(enz)*...
%                                     flux(strcmpi(model.rxns,'ACpts')); 
%     fout(it,tfr) = pvec.Vmax(tfr).*(1-1./(1+(pvec.KIact(fdp,tfr)./allmc(fdp)).^2));  
% 
% end
% M(1)  = 1;      % E
% M(2)  = 0.001;   % PEP
% M(3)  = 10;   % FBP

% 
% allhandles = feval(@Kotte2014glycolysis);
% rhsfunc = allhandles{2};
% func = @(t,x)rhsfunc(t,x,kEcat,KEacetate,...
%         KFbpFBP,vFbpmax,Lfbp,KFbpPEP,...
%         vEXmax,KEXPEP,...
%         vemax,KeFBP,ne,acetate,d);
% dMdt1 = func(0,M);  
% flux1 = Kotte_glycolysisflux(M,pvec);




plot(tout,yout);

plotKotteVariables(tout,yout,1);