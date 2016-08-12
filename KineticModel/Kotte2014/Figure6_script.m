% Figure 6
% change in kPEPout/bifurcation diagrams
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
    
idp = 11;
type = 'together';


% systems check
givenModel = @(t,x)KotteODE(t,x,model,pvec);
fluxg = Kotte_givenFlux([M;model.PM],pvec,model);
dMdtg = givenModel(0,M);

tspan = 0:0.1:6000;
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
    ap = 9;
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

%% plot bifurcation on acetate
% Figure 5 - load data from previous simulation
load('C:\Users\shyam\Documents\Courses\CHE1125Project\Results\KotteModel\kPEPoutVariation_Aug12.mat');

% needed variables: alliidpvec,alliidxeq,alliidfeq,tout,ap;
npts = size(alliidpvec,1);
nvar = size(alliidxeq,1);
ndp = size(alliidpvec,3);
tout = tspan;

for iid = 1:ndp
    if isfield(allnss,sprintf('iid%d',iid))
        msspts = find(allnss.(['iid' num2str(iid)]));
        sslps = allnss.(['iid' num2str(iid)])(msspts);
        ss = unique(sslps);
        nss = length(ss);
        allmsspts = [];
        for iss = 1:nss
            hf1 = [];
            allmsspts = union(allmsspts,msspts(sslps==ss(iss)));
            for ipt = 1:npts
                addanot.text = ['R' num2str(ipt)];
                % if point capable of mss
                if ismember(ipt,allmsspts)   
                    s1 =...
                    siid.(['iid' num2str(iid)]).(['pt' num2str(ipt)]).s1;
                    x1 =...
                    siid.(['iid' num2str(iid)]).(['pt' num2str(ipt)]).x1;
                    f1 =...
                    siid.(['iid' num2str(iid)]).(['pt' num2str(ipt)]).f1;
                    index =...
                    cat(1,siid.(['iid' num2str(iid)]).(['pt' num2str(ipt)]).s1.index);
                    hf1 = bifurcationPlot(x1,s1,f1,[4,1],ap,hf1,addanot);
                end
            end
        end
    end
end

% continue on kPEPout
% original system 
% get original ss and continuation without perturbations
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

% systems check
givenModel = @(t,x)KotteODE(t,x,model,pvec);
fluxg = Kotte_givenFlux([M;model.PM],pvec,model);
dMdtg = givenModel(0,M);    
    
opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
colorSpec = chooseColors(4,{'Green','Purple','Red','Orange'});
ac = find(strcmpi(model.mets,'ac[e]'));
npts = 1;

allxeq = zeros(length(M),npts);
allxdyn = zeros(length(M),length(tspan),npts);
allxf = zeros(length(M),npts);
allfeq = zeros(length(fluxg),npts);
allfdyn = zeros(length(fluxg),length(tspan),npts);
allpvec = pvec;

% conitnuation on kPEPout strating from 0
for ip = 1:length(idp)
    ap = 11;
    solveEquilibriumODE 
end

% continuation on kPEPout for perturbed values in cmb
hf2 = [];
ipt = 3;
while ipt<size(cmb,1)
% for ipt = 1:size(cmb,1)
    xeq = alliidxeq(:,ipt);
    pvec = alliidpvec(ipt,:);
    addanot.text = ['E' num2str(ipt)];
    for ip = 1:length(idp)
        ap = idp(ip);        
        % run MATCONT
        [data,y,p] = execMATCONT(xeq,pvec,ap,fluxg,model);
        if ~isempty(data) && size(data.s1,1)>2
            if ap == 11
                bifurcationPlot(data.x1,data.s1,data.f1,[4,1],ap);            
            end            
        end
        % save MATCONT results
        s.(['pt' num2str(ipt)]) = data;
    end
    ipt = ipt+4;
end