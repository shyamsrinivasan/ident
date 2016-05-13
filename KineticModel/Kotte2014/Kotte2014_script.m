% build stoichioemtrc matrices
addpath(genpath('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel'));
rxfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\Kotte2014\Kotte2014.txt';
cnfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\Kotte2014\Kotte2014C.txt';

% create model structure
[FBAmodel,parameter,variable,nrxn,nmetab] = modelgen(rxfname);

% obtain conentrations from file
[mc,FBAmodel,met] = readCNCfromFile(cnfname,FBAmodel);

% run FBA
Vup_struct.ACt2r = 1;
Vup_struct.ENZ1ex = 1;
FBAmodel = FBAfluxes(FBAmodel,'fba',{'ACt2r','ENZ1ex'},Vup_struct,...
                    [find(strcmpi(FBAmodel.rxns,'FDex'))...
                     find(strcmpi(FBAmodel.rxns,'PEPex'))]);
                 
% change bunds for FBAmodel and re run FBA
% [FBAmodel,bounds] = changebounds(FBAmodel,{'ACt2r','ENZ1ex'});
% FBAmodel.vl = bounds.vl;
% FBAmodel.vu = bounds.vu;
% 
% FBAmodel = FBAfluxes(FBAmodel,'fba',{'ACt2r','ENZ1ex'},Vup_struct,...
%                      find(strcmpi(FBAmodel.rxns,'EC_Biomass')));                 
                 
% flux envelope
[hfig,hsubfig,prxnid,flag] = FluxEnvelope(FBAmodel,...
                        {'bmt2r','PEPt2r';...
                        'ACpts','ENZtr';...
                        'EC_Biomass','ACpts';...
                        'ACpts','PEPt2r'},...
                        {'ACt2r','ENZ1ex'});
                    
% call to bifurcation analysis script using MATCONT
% KotteMATCONTscript         

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

% call to stoichioemtric analysis script
% Kotte_StoichiometricAnalysisScript

% call to parameter sampling script for analysis of mss
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
    
% sample parameters indicated by indices in idp
plb =  zeros(length(pvec),1);
pub = zeros(length(pvec),1);
idp = [1;4;7];
plb(idp) = [0.01;0.01;0.01];
pub(idp) = [100;100;100];
npts = 1000;
allpvec = sampleEKP(pvec,plb,pub,idp,npts);

% systems check
givenModel = @(t,x)KotteODE(t,x,model,pvec);
fluxg = Kotte_givenFlux([M;model.PM],pvec,model);
dMdtg = givenModel(0,M);

opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
tspan = 0:0.1:200;
[tout,yout] = ode45(givenModel,tspan,M,opts);
fout = zeros(length(tout),length(fluxg));
for it = 1:length(tout)
    fout(it,:) = Kotte_givenFlux([yout(it,:)';model.PM],pvec,model);
end
% plotKotteVariables(tout,yout,1);
% plotKotteVariables(tout,fout,2);

% run MATCONT on a multiple sets of parameters
% sample parameters indicated by indices in idp
% changing individual parameters
pvec = [kEcat,KEacetate,...
        KFbpFBP,vFbpmax,Lfbp,KFbpPEP,...
        vEXmax,KEXPEP,...
        vemax,KeFBP,ne,acetate,d];
idp = [1;4;7];
vals = [0.01 0.1 1 10 100];
npts = 500;
allpvec = discreteEKPsample(pvec,vals,idp,npts);

xeq = yout(end,:)';
runMATCONT % run MATCONT from script
sinit = s1;
xinit = x1;
finit = f1;

% get the mss for y and p
[yss,iyval,fyval] = parseMATCONTresult(sinit,y);
[pss,ipval,fpval] = parseMATCONTresult(sinit,p);

nss = size(yss,2);
flux1 = zeros(5,nss);
for iss = 1:nss
    pvec(ap) = p(iss);
    flux1(:,iss) = KotteMATCONTflux(yss(:,iss),pvec);
end

plotPointsonFluxEnvelope(hfig,hsubfig,[3 5;1 2;3 1;1 5],flux1)

% use CNA to calculate EFMs and plot them on the envelope
% call to stoichioemtric analysis script
Kotte_StoichiometricAnalysisScript

% Caluclate yield with MCS
nmcs = size(mcs,1);


% use ADMAT to calculate jacobians
admatfun = @(x)Kotte_givenNLAE(x,model,pvec);
% x = ones(length(M),1);
xADMATobj = deriv(M,eye(3));
xADMATres = admatfun(xADMATobj);
F = getval(xADMATres);
J = getydot(xADMATres);


Jxact = KottegivenJacobian(M,pvec,model);

% permutations
ps = nchoosek([0.01 0.1 1 10 100],3);
vs = zeros(3,0);
for ip = 1:size(ps,1)
    vp = perms(ps(ip,:));
    vs = [vs vp'];
end
vs = unique(vs','rows');
vs = [0.01 0.01 0.01;.1 .1 .1;10 10 10;100 100 100];

npts = size(vs,1);

tspan = 0:0.1:50000;
allyoutss = zeros(length(M),npts);
allxeq = zeros(length(M),npts);
allxf = zeros(length(M),npts);
allflag = zeros(1,npts);
for ipt = 1:npts
    fprintf('Iteration #%d Equilibrium Integration...',ipt);
    % change in pvec
%     pvec = allpvec(ipt,:);
    pvec(idp) = vs(ipt,:);
    
    % new equilibrium solution
    givenModel = @(t,x)KotteODE(t,x,model,pvec);
    [tout,yout] = ode45(givenModel,tspan,M,opts);
    allyoutss(:,ipt) = yout(end,:)';    
    plotKotteVariables(tout,yout,1);    
    
    gfun = @(x)Kotte_givenNLAE(x,model,pvec);
    options = optimoptions('fsolve','Display','iter',...
                                    'TolFun',1e-12,'TolX',1e-12,...
                                    'MaxFunEvals',10000,...
                                    'MaxIter',5000);
%     [xf,fval,exitflag,output,jacobian] = fsolve(gfun,M,options);
    xeq = allyoutss(:,ipt);
    allxeq(:,ipt) = xeq;
%     allxf(:,ipt) = xf;
%     allflag(1,ipt) = exitflag;
%     xeq = xf;
    
    fprintf('Complete\n');
    
    % continuation from initial equilibrium - initialization
    fprintf('Iteration #%d Equilibrium Continuation...',ipt);
    % run MATCONT
    runMATCONT
    
    % save MATCONT results
    s.(['pt' num2str(ipt)]).s1 = s1;
    s.(['pt' num2str(ipt)]).x1 = x1;
    s.(['pt' num2str(ipt)]).f1 = f1;
    
    fprintf('Complete\n');
end

% check which solutions have mss
mssid = [];
for ipt = 1:npts
    s1 = s.(['pt' num2str(ipt)]).s1;
    nLP = size(s1,1);
    if nLP > 2
        fprintf('Vector %d has %d Steady States\n',ipt,nLP);
        mssid = union(mssid,ipt);
    end
end

% calculation of fluxes for allxeq
allfeq = zeros(length(fluxg),npts);
for ipt = 1:npts
    pvec(idp) = vs(ipt,:);
    allfeq(:,ipt) = Kotte_givenFlux([allxeq(:,ipt);model.PM],pvec,model);
end

% plot solutions that have mss


% ode for different parameter sets
% vary all parameters simulataneously
% allyout = zeros(length(tspan),length(M),npts);
% allyoutss = zeros(length(M),npts);
% allfout = zeros(length(tspan),length(fluxg),npts);
% for ipt = 1:npts
%     fprintf('Iteration #%d...',ipt);
%     pvec = allpvec(ipt,:);
%     givenModel = @(t,x)KotteODE(t,x,model,pvec);
%     [tout,yout] = ode45(givenModel,tspan,M,opts);
%     allyout(:,:,ipt) = yout;
%     allyoutss(:,ipt) = yout(end,:)';
%     fprintf('Complete\n');
% end
% plotKotteVariables(tout,allyout,1);
% plotKotteVariables(allpvec(:,idp)',allyoutss,3);

  
    
    

% mfn = CNAloadNetwork(1,true,true);

% FBAmodel.rxns(cellfun(@(x)strcmpi(x,'EC_Biomass'),FBAmodel.rxns)) = {'mue'};






% ffca - feasibility flux coupling analysis
% network.stoichiometricMatrix = FBAmodel.S;
% network.reversibilityVector = FBAmodel.rev;
% network.Reactions = FBAmodel.rxns;
% network.Metabolites = FBAmodel.mets;
% 
% [fctable,blocked] = FFCA('glpk',network);

% efms using METATOOLS
% rd_ems = nsa_em(full(FBAmodel.S),~FBAmodel.rev');



% Kotte_Cscript
% allhandles = feval(@Kotte2014Ckinetics);
% rhsfunc = allhandles{2};
% CModel = @(t,x)rhsfunc(t,x,model,kEcat,KEacetate,...
%         KFbpFBP,vFbpmax,Lfbp,KFbpPEP,...
%         vEXmax,KEXPEP,...
%         vemax,KeFBP,ne,acetate,d);
% fluxC = Kotte_CFlux([x1;model.PM],pvec,model);
% dMdtC = CModel(0,x1);

% C model fsolve
% cfun = @(x)Kotte_CNLAE(x,model,pvec);
% dMc = cfun(x1);
% options = optimoptions('fsolve','Display','iter',...
%                        'TolFun',1e-10,...
%                        'TolX',1e-10,...
%                        'MaxFunEvals',1000000,...
%                        'MaxIter',50000);
% [x2,fval,exitflag,output,jacobian] = fsolve(cfun,x1,options);
% fcout = Kotte_CFlux([M;model.PM],pvec,model);    
% 
% % C model SS
% opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
% [tout,yout] = ode45(CModel,0:0.1:60,M,opts);
% fout = zeros(length(tout),4);
% for it = 1:length(tout)
%     fout(it,:) = Kotte_givenFlux([yout(it,:)';model.PM],pvec,model);
% end

% pvec.d = 0.25;
% allhandles = feval(@Kotte2014Ckinetics);
% rhsfunc = allhandles{2};
% func = @(t,x)KotteCkinetics(t,x,model,pvec);
% dMdt = func(0,M);
% allmc = [M;model.PM];
% 
% fun = @(x)KotteCkinetics(x,pvec,model);
% dMdt = fun(M);
% 
% 
% [x1,fval,exitflag,output,jacobian] = fsolve(fun,M,options);

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




