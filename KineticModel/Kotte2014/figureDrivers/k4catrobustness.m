% sampling from edge of bistable region
% the edge of the region of bistability is the most desirable point for
% systems to be at in terms of production capabilities
% what happens to systems sampled around this point?
% how robust are systems sampled around the this point to reaching either
% the high state or the low state

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
[model,~,newmc,cnstmet] =...
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
k1cat = 1;
K1ac = 0.1;    % or 0.02
K3fdp = 0.1;
v3max = 1;
L3 = 4e6;
K3pep = 0.1;
v2max = 1;
K2pep = 0.3;
vemax = 1.1;        % for bifurcation analysis: 0.7:0.1:1.3
KeFBP = 0.45;       % or 0.45
ne = 2;             % or 2
acetate = 0.1;      % a.u acetate
d = 0.25;           % or 0.25 or 0.35
k4cat = 0.2;
pvec = [K1ac,K3fdp,L3,K3pep,...
        K2pep,vemax,KeFBP,ne,acetate,d,...
        k4cat,k1cat,v3max,v2max];

% systems check
givenModel = @(t,x)KotteODE(t,x,model,pvec);
fluxg = Kotte_givenFlux([M;model.PM],pvec,model);
tspan = 0:0.1:2000;
opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
ac = find(strcmpi(model.mets,'ac[e]'));
npts = 1;
ap = 9;
allpvec = pvec;

% find equilibrium solution and run equilibrium continuation
allxeq = zeros(length(M),npts);
allxdyn = zeros(length(M),length(tspan),npts);
allxf = zeros(length(M),npts);
allfeq = zeros(length(fluxg),npts);
allfdyn = zeros(length(fluxg),length(tspan),npts);
solveEquilibriumODE

% get saddle node to get both stable steady states and get a bistable
% acetate concentration for perturbation
[orig_saddle,orig_saddlepar] = getsaddlenode(data.s1,data.x1,5e-3);
pvec(ap) = orig_saddlepar;
model.PM(ac-length(orig_saddle)) = orig_saddlepar;

% perturb saddle to get steady states
eps = 1e-4;
tspanf = 0:0.1:2000;
pival = orig_saddle+eps*[1;1;1];
[~,xeq1,~,feq1] = solveODEonly(1,pival,model,pvec,opts,tspanf);
nival = orig_saddle-eps*[1;1;1];
[~,xeq2,~,feq2] = solveODEonly(1,nival,model,pvec,opts,tspanf);
xss = [xeq1 xeq2];

idp = 11;
type = 'together';
cmb = linspace(0,1,50)';
npts = size(cmb,1);

% set acetate conentration to saddle
pvec = [K1ac,K3fdp,L3,K3pep,...
        K2pep,vemax,KeFBP,ne,acetate,d,...
        k4cat,k1cat,v3max,v2max];
% pvec(ap) = orig_saddlepar;
% model.PM(ac-length(orig_saddle)) = orig_saddlepar;

% set parameters from cmb at idp position(s)
allpvec = repmat(pvec,npts,1);
allpvec(:,idp) = cmb;

% find equilibrium point for lowest acetate 
allpvec(:,ap) = 0.01;
model.PM(ac-length(orig_saddle)) = 0.01;    
[~,allxeqlac] = solveODEonly(npts,M,model,allpvec,opts,tspan);

% and continue on acetate
[s,mssid,nss] = setupMATCONT(allxeqlac,allpvec,ap,model,fluxg,npts,1500);

%% get pep and v4 vs k4cat for different acetate 
acetate = orig_saddlepar;
ap = 9;
colorSpec = chooseColors(5,{'Green','Purple','Red','Navy','HotPink'});
saddleac = zeros(npts,length(acetate));
xeqac = zeros(2*nvar,npts,length(acetate));
feqac = zeros(2*length(fluxg),npts,length(acetate));
for iac = 1:length(acetate)    
    % calculate saddle for each acetate concentration
    eps = 1e-4;
    [saddle,saddlepar,status] = eqptwrapper(s,nvar,acetate(iac),eps);
    
    % good saddle node points only
    goodsaddle = saddle(:,logical(status));
    goodsaddlepar = saddlepar(logical(status));
    saddleac(logical(status),iac) = goodsaddlepar;
    allpvec(logical(status),ap) = goodsaddlepar;
    % saddle parameter out of bifurcation bounds
    % get the one possible steady state
    oubsaddle = saddle(:,~logical(status));
    oubsaddlepar = saddlepar(~logical(status));
    saddleac(~logical(status),iac) = acetate(iac);
    allpvec(~logical(status),ap) = acetate(iac);
    
    for ipt = 1:npts
        if ismember(ipt,find(status))
            model.PM(ac-length(orig_saddle)) = saddlepar(ipt); 
            % perturb saddle to get steady states
            eps = 1e-4;                            
            pival = saddle(:,ipt)+eps*[1;1;1];
            [~,xeq1,~,feq1] =...
            solveODEonly(1,pival,model,allpvec(ipt,:),opts,tspanf);
            nival = saddle(:,ipt)-eps*[1;1;1];
            [~,xeq2,~,feq2] =...
            solveODEonly(1,nival,model,allpvec(ipt,:),opts,tspanf);            
            xeqac(nvar+1:end,ipt,iac) = xeq2;
            feqac(length(fluxg)+1:end,ipt,iac) = feq2;
        else
            % get the only possible steady state
            model.PM(ac-length(orig_saddle)) = acetate(iac);
            [~,xeq1,~,feq1] = solveODEonly(1,M,model,allpvec(ipt,:),opts,tspan);            
            xeqac(nvar+1:end,ipt,iac) = xeq1;
        end
        xeqac(1:nvar,ipt,iac) = xeq1;
        feqac(1:length(fluxg),ipt,iac) = feq1;
    end
    
    xeqid = zeros(2,npts);
    xeqid(1,xeqac(1,:)>xeqac(2,:)) = 1;
    xeqid(1,xeqac(2,:)>xeqac(1,:)) = 2;
    xeqid(2,xeqac(4,:)>xeqac(5,:)) = 1;
    xeqid(2,xeqac(5,:)>xeqac(4,:)) = 2;    
    % get indices where xeqid(1) != xeqid(2)
    shiftid = find(xeqid(1,:)~=xeqid(2,:),1,'last');
    % adjust pvec?
    pvec = allpvec(shiftid,:);
    model.PM(ac-nvar) = pvec(9);
    % sample points around shiftid (ball of radius?)
    samples = SphereSampling(xeqac(1:3,shiftid),1,1000);
    samples = samples(:,samples(1,:)>0);
    samples = samples(:,samples(2,:)>0);
    samples = samples(:,samples(3,:)>0);
    % simulate from each sample
    nsample = size(samples,2);
    smpxeq = zeros(nvar,nsample);
    smpfeq = zeros(length(fluxg),nsample);
    %%
    parfor isample = 1:nsample
        [~,xeqs,~,feqs] =...
        solveODEonly(1,samples(:,isample),model,pvec,opts,tspanf);
        smpxeq(:,isample) = xeqs;
        smpfeq(:,isample) = feqs;
    end    
    % find edge of bistable region for each acetate 
    % point at which status changes a second time?
    % point at which you stop get 2 steady states a second time?
    % find points between changes in status
end

