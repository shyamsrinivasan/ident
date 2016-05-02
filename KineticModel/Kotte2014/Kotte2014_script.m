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
                 
% flux envelope
[hsubfig,prxnid,flag] = FluxEnvelope(FBAmodel,...
                        {'bmex','PEPex'},...
                        {'ACt2r','ENZ1ex'});

% call to bifurcation analysis script using MATCONT
% KotteMATCONTscript
                    
% remove metabolites held constant from consideration in the model
% integration phase
[model,pvec,newmc,cnstmet] =...
remove_eMets(FBAmodel,parameter,mc,[FBAmodel.Vind FBAmodel.Vex],...
{'enz1[c]','enz1[e]','enz[e]','ac[e]','bm[c]','bm[e]','pep[e]'});

% change bunds for FBAmodel
[FBAmodel,bounds] = changebounds(FBAmodel,{'ACt2r','ENZ1ex'});
FBAmodel.vl = bounds.vl;
FBAmodel.vu = bounds.vu;

FBAmodel = FBAfluxes(FBAmodel,'fba',{'ACt2r','ENZ1ex'},Vup_struct,...
                     find(strcmpi(FBAmodel.rxns,'EC_Biomass')));

% EM analysis using Cell Net Analyzer (CNA)
% convert model to conform to CNA form

% mfn = CNAloadNetwork(1,true,true);

% FBAmodel.rxns(cellfun(@(x)strcmpi(x,'EC_Biomass'),FBAmodel.rxns)) = {'mue'};
spec = ones(size(FBAmodel.S,1),1)';
spec(1:FBAmodel.nint_metab) = 0;
cnap.has_gui = 0;
cnap.net_var_name = 'KotteGmodel';
cnap.type = 1;
cnap.specID = char(FBAmodel.mets); 
cnap.specLongName = char(FBAmodel.mets);
cnap.specExternal = spec;
cnap.specInternal = find(~cnap.specExternal);
cnap.nums = size(FBAmodel.S,1);
cnap.numis = size(cnap.specInternal,2);
% cnap.macroID = 'BC1';
% cnap.macroLongName = 'BC1';
% cnap.macroComposition =...
% sparse(find(strcmpi(FBAmodel.mets,'bm[c]')),1,1,length(FBAmodel.mets),1);
% cnap.macroDefault = 1;
% cnap.nummac = 1;
cnap.stoichMat = full(FBAmodel.S); % zeros(length(FBAmodel.mets),1)];
cnap.numr = size(cnap.stoichMat,2);
cnap.reacID = char(FBAmodel.rxns); % char([FBAmodel.rxns;'mue']);
cnap.objFunc = zeros(cnap.numr,1);
cnap.reacMin = FBAmodel.vl; % 0];
cnap.reacMax = FBAmodel.vu; % 100];


[cnap,errval] = CNAgenerateMFNetwork(cnap);

cnap.path =...
'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\Kotte2014\KotteGmodel';
cnap = CNAsaveNetwork(cnap);

% flux optimization
% constr = zeros(cnap.numr,1);
% constr(constr==0) = NaN;
% constr(1) = 10;
% constr(4) = 10;
% constr(9) = 0;
% constr(14) = 1;
% [flux,success,status] = CNAoptimizeFlux(cnap,constr,cnap.macroDefault,2,2);

% flux variablity
% reacval = zeros(cnap.numr,1);
% reacval(reacval==0) = NaN;
% reacval(5) = -10;
% reacval(11) = -10;
% [minFlux,maxFlux,success,status] =...
% CNAfluxVariability(cnap,reacval,cnap.macroDefault,2);

% remove conserved quantities
% [cnap,delspec] = CNAremoveConsRel(cnap,1,0,0);

% phase plane analysis
% status = CNAplotPhasePlane(cnap,reacval,cnap.macroDefault,[14;10;12;13],2);

% EFM calculation
constr = zeros(cnap.numr,4);
constr(constr==0) = NaN;
% exclude efms w/ exchanges only
constr([5 10 11 12 13],1) = 0;
% include efms with only all major reactions
% constr([8],1) = 1;
% constr(1,2) = 1;
% constr(9,2) = 0.1961;
constr(14,2) = 0.3;
constr(8,2) = 1;
constr(1,3) = 1;
% constr(7,2) = -1;
% % constr(14,2) = 0.1;
% constr(9,2) = 0.1;
% constr([5 7],3) = 0;
% constr([1 4],4) = 1;

% [efm,rev,idx,ray] = CNAcomputeEFM(cnap,constr,2,1,0,0);
% printEFM(efm,idx,ray,cnap);

% reacID = cnap.reacID(idx,:);

% cut set calculation
% target = efm([1 4],:);
% set2save(1).tabl2save = efm([1 2],:);
% set2save(2).tabl2save = efm(2,:);
% set2save(1).min2save = 2;
% set2save(2).min2save = 1;
% set2save = [];
% cutsets = CNAcomputeCutsets(target,10,reacID,set2save);
% printCS(cutsets,reacID);

% constrained minimal cut sets - 
% ends with an  error due to lack of MCSs
cnap.reacMin(cnap.reacMin == -100) = -Inf;
cnap.reacMax(cnap.reacMax == 100) = Inf;

tar = zeros(1,cnap.numr);
nT = 1;
T = repmat(tar,nT,1);
% 
cellreacID = cellstr(cnap.reacID);
cellreacID = cellfun(@(x)strtrim(x),cellreacID,'UniformOutput',false);
T(1,strcmpi(cellreacID,'PEPt2r')) = -1;
% T(3,strcmpi(cellreacID,'ACpts')) = -1;
t = zeros(nT,1);
% t(1) = 0.5;
t(1) = -0.2;
% t(3) = -1;

nD = 1;
D = repmat(tar,nD,1);
D(1,strcmpi(cellreacID,'ACpts')) = -1;
% D(2,strcmpi(cellreacID,'PEPt2r')) = -1;
% D(3,strcmpi(cellreacID,'ACpts')) = 1;
d = zeros(nD,1);
d(1) = 0;
% d(2) = 0;
% d(3) = 1;

notknockable = [];
maxMCS = 100;
maxMCSsize = 5;
filename = [];

mcs =...
CNAMCSEnumerator(cnap,T,t,[],[],notknockable,maxMCS,maxMCSsize,filename);
printCS(mcs,cnap.reacID);
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




