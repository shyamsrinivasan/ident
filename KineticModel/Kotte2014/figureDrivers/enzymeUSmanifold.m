% get new unstable manifolds for different enzyme expression levels
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

% run enzyme perturbations
% sample parameters indicated by indices in idp
cmb = [.05 1 1;1 .05 1;1 1 .05;.05 .05 .05;...
       .125 1 1;1 .125 1;1 1 .125;.125 .125 .125;...
       .25 1 1;1 .25 1;1 1 .25;.25 .25 .25;...
       .5 1 1;1 .5 1;1 1 .5;.5 .5 .5;...
       2 1 1;1 2 1;1 1 2;2 2 2;...
       4 1 1;1 4 1;1 1 4;4 4 4];

idp = [12 13 14];
type = 'together';
npts = size(cmb,1);

if strcmpi(type,'together')
    alliidpvec = zeros(npts,length(pvec),size(idp,1));    
else    
    alliidpvec = zeros(npts,length(pvec),length(idp));    
end

% set acetate conentration
pvec = [K1ac,K3fdp,L3,K3pep,...
        K2pep,vemax,KeFBP,ne,acetate,d,...
        k4cat,k1cat,v3max,v2max];
pvec(ap) = 0.01;
model.PM(ac-length(orig_saddle)) = 0.01;

% set parameters from cmb at idp position(s)
allpvec = repmat(pvec,npts,1);
allpvec(:,idp) = cmb;

% find equilibrium point for lowest acetate 
[~,allxeqlac] = solveODEonly(npts,M,model,allpvec,opts,tspan);

% and continue on acetate
[s,mssid,nss] = setupMATCONT(allxeqlac,allpvec,ap,model,fluxg,npts,1500);

%% get boundaries of acetate bistability
allEbsval = allpvec(mssid,idp);
acbounds = zeros(2,length(mssid)); % [min;max];
% xbounds = zeros(nvar,2*length(mssid));
mssipt = 1;
for ipt = 1:npts
    if ismember(ipt,mssid)
        index = cat(1,s.(['pt' num2str(ipt)]).s1.index);
        x1 = s.(['pt' num2str(ipt)]).x1;
%         xcont = x1(1:nvar,index);
        pcont = x1(nvar+1:end,index);
%         xbounds(:,2*mssipt-1:2*mssipt) = xcont(:,2:end-1);
        acbounds(1,mssipt) = min(pcont(:,2:end-1));
        acbounds(2,mssipt) = max(pcont(:,2:end-1));
        mssipt = mssipt+1;
    end
end

%% get a saddle node for each of the bistable cases and plot Wus
alleig = zeros(nvar,length(mssid));
allw = zeros(nvar,nvar,length(mssid));
allsaddles = zeros(nvar,length(mssid));
allxeq = zeros(nvar,length(mssid));
issid = 1;
tspanr = 0:-.1:-30;
tspanf = 0:0.1:2000;
eps = 1e-4;
eps1 = 1e-4;
allWus = zeros(nvar,length(tspanr),length(mssid));
% pick a saddle node    
[saddle,saddlepar,status] = geteqcontpt(s,orig_saddlepar,eps);

figure
hold on
for iss = 1:npts
    if ismember(iss,mssid)
        ispvec = allpvec(iss,:);
                
        % get eig val and eig vector
        ispvec(ap) = saddlepar(iss);
        model.PM(ac-nvar) = saddlepar(iss);
        [~,alleig(:,issid),w] = getKotteJacobian(saddle(:,iss),ispvec,model);
        allw(:,:,issid) = w;
        
        % get Wus at saddlepar        
        [xWus,xeq] =...
        calc1DWus(saddle(:,iss),w,alleig(:,issid),model,ispvec,opts,tspanr,tspanf,eps1);
        
        % collect reverse trajectory 
        allWus(:,:,issid) = xWus;
        allxeq(:,issid) = xeq(:,1);      
        
        issid = issid+1;
    end
end

% unstable manifolds after chopping
relallWus = real(allWus);
figure
hold on
for jpt = 1:length(mssid)
    [~,nzid,~] = find(relallWus(:,:,jpt)~=0,1,'last');
    relWus = relallWus(:,1:nzid,jpt);
    kid = relWus(1,:)<0|relWus(1,:)>20|relWus(2,:)<0|relWus(2,:)>20|relWus(3,:)<0|relWus(3,:)>20;
    relWus(:,kid) = [];
%     allx = [allx relWus(1,:)];
%     ally = [ally relWus(2,:)];
%     allz = [allz relWus(3,:)];
    % plot trajectory
    plot3(relWus(1,:),relWus(2,:),relWus(3,:),'LineWidth',2);
    plot3(saddle(1,mssid(jpt)),saddle(2,mssid(jpt)),saddle(3,mssid(jpt)),'Marker','.','Color','r','MarkerSize',16);
    drawnow
end

%% extract all unstable points for each enzyme perturbation that is bistable
% figure
hold on
for ipt = 1:npts
    if ismember(ipt,mssid)
        allsaddles = extractsaddles(s.(['pt' num2str(ipt)]).s1,s.(['pt' num2str(ipt)]).x1);       
        nsadpts = size(allsaddles,2);
        plot3(allsaddles(1,:),allsaddles(2,:),allsaddles(3,:),'Color','r');
    end
end



%%
colorSpec = chooseColors(4,{'Navy','HotPink','Red','Orange'});
for iid = 1:1 % length(idp)
    fprintf('Parameter Combination #%d\n',iid); 
    % find equilibrium solution from different ss at the original saddle
    % node
    alliidpvec(:,:,iid) = allpvec;    
    hf1 = [];
    ha1 = [];
    hf2 = [];
    ha2 = [];

    ivalpts = zeros(2*nvar,npts);
    xeqpts = zeros(2*nvar,npts);
    eqid = zeros(2,npts);
    ivalid = zeros(2,npts);
    % perturbation for all points
    for ipt = 1:npts
        pvec = alliidpvec(ipt,:,iid);                
        % perturbations from ss (xeq1 and xeq2)
        [ivalpts,ivalid,xeqpts,eqid,hf1,ha1] = ParameterPerturbations(model,pvec,...
            xss,ivalpts,ivalid,xeqptseqid,ipt,tspanf,colorSpec,opts,hf1,ha1);                          
    end    
end