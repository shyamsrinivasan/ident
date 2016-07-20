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
% [hfig,hsubfig,prxnid,flag] = FluxEnvelope(FBAmodel,...
%                         {'bmt2r','PEPt2r';...
%                         'ACpts','ENZtr';...
%                         'bmt2r','ACpts';...
%                         'ACpts','PEPt2r'},...
%                         {'ACt2r','ENZ1ex'});
                    
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

% use CNA to calculate EFMs and plot them on the envelope
% change model to be CNA compliant
% [cnap,errval] = getCNAmodel(FBAmodel);
% call to stoichioemtric analysis script for calculating EFMs
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
kPEPout = 0.2;
pvec = [KEacetate,KFbpFBP,Lfbp,KFbpPEP,...
        KEXPEP,vemax,KeFBP,ne,acetate,d,...
        kPEPout,kEcat,vFbpmax,vEXmax];

% npts = 1;
% allpvec = pvec;
% solveEquilibriumODE
    
% permutations of parameter values 
% ps = nchoosek([0.01 0.1 1 10 100],3);
% vs = zeros(3,0);
% for ip = 1:size(ps,1)
%     vp = perms(ps(ip,:));
%     vs = [vs vp'];
% end
% vs = unique(vs','rows');
% vs = [0.001 0.001 0.001;0.01 0.01 0.01;.1 .1 .1;1 1 1;10 10 10];
% npts = size(vs,1);

% sample parameters indicated by indices in idp
cmb = [1e6;2e6;3e6;4e6;5e6;6e6;7e6;8e6;9e6;1e7];
% cmb = [.05 1 1;1 .05 1;1 1 .05;.05 .05 .05;...
%        .25 1 1;1 .25 1;1 1 .25;.25 .25 .25;...
%        .5 1 1;1 .5 1;1 1 .5;.5 .5 .5;...
%        2 1 1;1 2 1;1 1 2;2 2 2;...
%        4 1 1;1 4 1;1 1 4;4 4 4];
idp = 3;
type = 'together';
npts = size(cmb,1);

% systems check
givenModel = @(t,x)KotteODE(t,x,model,pvec);
fluxg = Kotte_givenFlux([M;model.PM],pvec,model);
dMdtg = givenModel(0,M);

tspan = 0:0.1:2000;
% allflag = zeros(1,npts);

opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
if strcmpi(type,'together')
    alliidpvec = zeros(npts,length(pvec),size(idp,1));
    alliidxeq = zeros(length(M),npts,size(idp,1));
    alliidxdyn = zeros(length(M),length(tspan),npts,size(idp,1));
    alliidfeq = zeros(length(fluxg),npts,size(idp,1));
    alliidfdyn = zeros(length(fluxg),length(tspan),npts,size(idp,1));
else    
    alliidpvec = zeros(npts,length(pvec),length(idp));
    alliidxeq = zeros(length(M),npts,length(idp));
    alliidxdyn = zeros(length(M),length(tspan),npts,length(idp));
    alliidfeq = zeros(length(fluxg),npts,length(idp));
    alliidfdyn = zeros(length(fluxg),length(tspan),npts,length(idp));
end
allpvec = repmat(pvec,npts,1);
allpvec(:,idp) = cmb;
   
for iid = 1:1 % length(idp)
    % reset pvec
    pvec = [KEacetate,KFbpFBP,Lfbp,KFbpPEP,...
            KEXPEP,vemax,KeFBP,ne,acetate,d,...
            kPEPout,kEcat,vFbpmax,vEXmax];
    
    plb = 0;
    pub = 1;    
    fprintf('Parameter Combination #%d\n',iid);
    
%     allpvec = sampleEKP(pvec,plb,pub,idp(iid),npts); 
%     smp_pvec = linspace(0,1,5000);
%     allpvec(:,idp(iid)) = smp_pvec';
    
    % run equilibrium solution followed by MATCONT
    allxeq = zeros(length(M),npts);
    allxdyn = zeros(length(M),length(tspan),npts);    
    allfeq = zeros(length(fluxg),npts);
    allfdyn = zeros(length(fluxg),length(tspan),npts);
    solveEquilibriumODE
    
    % save solution
    alliidpvec(:,:,iid) = allpvec;
    alliidxeq(:,:,iid) = allxeq;
    alliidxdyn(:,:,:,iid) = allxdyn;
    alliidfeq(:,:,iid) = allfeq;
    alliidfdyn(:,:,:,iid) = allfdyn;
    
    siid.(['iid' num2str(iid)]) = s;
    allmssid.(['iid' num2str(iid)]) = mssid;
    allnss.(['iid' num2str(iid)]) = nss;
end


% tspan = 0:0.1:200;
% [tout,yout] = ode45(givenModel,tspan,M,opts);
% fout = zeros(length(tout),length(fluxg));
% for it = 1:length(tout)
%     fout(it,:) = Kotte_givenFlux([yout(it,:)';model.PM],pvec,model);
% end
% plotKotteVariables(tout,yout,1);
% plotKotteVariables(tout,fout,2);

% use CNA to calculate EFMs and plot them on the envelope
% change model to be CNA compliant
% [cnap,errval] = getCNAmodel(FBAmodel);
% call to stoichioemtric analysis script for calculating EFMs
% Kotte_StoichiometricAnalysisScript

% Caluclate yield with MCS
% nmcs = size(mcs,1);


% use ADMAT to calculate jacobians
% admatfun = @(x)Kotte_givenNLAE(x,model,pvec);
% % x = ones(length(M),1);
% xADMATobj = deriv(M,eye(3));
% xADMATres = admatfun(xADMATobj);
% F = getval(xADMATres);
% J = getydot(xADMATres); 
% Jxact = KottegivenJacobian(M,pvec,model);




% convert Kotte flux to model flux

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




