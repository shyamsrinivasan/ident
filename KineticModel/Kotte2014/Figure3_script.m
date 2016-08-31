% Figure 3
% calculate changes in steady state due to parameter pertrubations
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
    
% systems check
fluxg = Kotte_givenFlux([M;model.PM],pvec,model);
tspan = 0:0.1:2000;
opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
ac = find(strcmpi(model.mets,'ac[e]'));
npts = 1;
ap = 9;

allxeq = zeros(length(M),npts);
allxdyn = zeros(length(M),length(tspan),npts);
allxf = zeros(length(M),npts);
allfeq = zeros(length(fluxg),npts);
allfdyn = zeros(length(fluxg),length(tspan),npts);
allpvec = pvec;

% find equilibrium solution and run equilibrium continuation
solveEquilibriumODE;          

% get saddle node
[orig_saddle,orig_saddlepar] = getsaddlenode(data.s1,data.x1,5e-3);
%% run enzyme perturbations
% sample parameters indicated by indices in idp
cmb = [.05 1 1;1 .05 1;1 1 .05;.05 .05 .05;...
       .125 1 1;1 .125 1;1 1 .125;.125 .125 .125;...
       .25 1 1;1 .25 1;1 1 .25;.25 .25 .25;...
       .5 1 1;1 .5 1;1 1 .5;.5 .5 .5;...
       2 1 1;1 2 1;1 1 2;2 2 2;...
       4 1 1;1 4 1;1 1 4;4 4 4];
% e_exp = linspace(0.005,4,50);
% cmb = ones(50*4,3);
% i = 0;
% while i<50
%     cmb(1+4*i,1) = e_exp(i+1);
%     cmb(2+4*i,2) = e_exp(i+1);
%     cmb(3+4*i,3) = e_exp(i+1);
%     cmb(4+4*i,1:3) = repmat(e_exp(i+1),1,3);
%     i = i+1;
% end   
idp = [12 13 14];
type = 'together';
npts = size(cmb,1);

% systems check
givenModel = @(t,x)KotteODE(t,x,model,pvec);
fluxg = Kotte_givenFlux([M;model.PM],pvec,model);

tspan = 0:0.1:2000;

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

% set acetate conentration
pvec = [KEacetate,KFbpFBP,Lfbp,KFbpPEP,...
        KEXPEP,vemax,KeFBP,ne,acetate,d,...
        kPEPout,kEcat,vFbpmax,vEXmax];
pvec(ap) = orig_saddlepar;
model.PM(ac-length(orig_saddle)) = orig_saddlepar;    

% set parameters from cmb at idp position(s)
allpvec = repmat(pvec,npts,1);
allpvec(:,idp) = cmb;

% save allpvec for ap
allpvecofap = allpvec(:,ap);

for iid = 1:1 % length(idp)
    fprintf('Parameter Combination #%d\n',iid);    
    
    % find equilibrium solution followed by MATCONT
    allxeq = zeros(length(M),npts);
    allxdyn = zeros(length(M),length(tspan),npts);    
    allfeq = zeros(length(fluxg),npts);
    allfdyn = zeros(length(fluxg),length(tspan),npts);
    [allxdyn,allxeq,allfdyn,allfeq] =...
    solveODEonly(npts,M,model,allpvec,opts,tspan,...
              allxdyn,allxeq,allfdyn,allfeq);
          
    % continue on acetate for all equilibirum solutions to different
    % parameter combinations
    ap = 9;
%     allpvec(:,ap) = 0.01;
    [s,mssid,nss] = setupMATCONT(allxeq,allpvec,ap,model,fluxg,npts);
%     solveEquilibriumODE
    
    % restore allpvec(ap) after continuation
    allpvec(:,ap) = allpvecofap;
    
    % save solution
    alliidpvec(:,:,iid) = allpvec;
    alliidxeq(:,:,iid) = allxeq;
    alliidxdyn(:,:,:,iid) = allxdyn;
    alliidfeq(:,:,iid) = allfeq;
    alliidfdyn(:,:,:,iid) = allfdyn;
    
    siid.(['iid' num2str(iid)]) = s;
    allmssid.(['iid' num2str(iid)]) = mssid;
    allnss.(['iid' num2str(iid)]) = nss;
    
    % reset pvec for next iteration of iid - need to check before using
    pvec = [KEacetate,KFbpFBP,Lfbp,KFbpPEP,...
            KEXPEP,vemax,KeFBP,ne,acetate,d,...
            kPEPout,kEcat,vFbpmax,vEXmax];
    pvec(ap) = orig_saddlepar;    
end

%% runDynamicRep
% Figure 3
load('C:\Users\shyam\Documents\Courses\CHE1125Project\Results\KotteModel\VmaxVariation_EquilibriumData_Aug30.mat');

% get figure for para,eter perturbation through initial value perturbation
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

opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
colorSpec = chooseColors(4,{'Navy','HotPink','Red','Orange'});
ac = find(strcmpi(model.mets,'ac[e]'));
npts = 1;
ap = 9;

allxeq = zeros(length(M),npts);
allxdyn = zeros(length(M),length(tspan),npts);
allxf = zeros(length(M),npts);
allfeq = zeros(length(fluxg),npts);
allfdyn = zeros(length(fluxg),length(tspan),npts);
allpvec = pvec;
solveEquilibriumODE     

% get saddle node
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
fss = [feq1 feq2];
xss = [xeq1 xeq2];

% needed variables: alliidpvec,alliidxeq,alliidfeq,tout,ap;
npts = size(alliidpvec,1);
nvar = size(alliidxeq,1);
ndp = size(alliidpvec,3);
tout = tspan;

% pertrubation calculation for all paramtere combinaions
for iid = 1:ndp
    hf1 = [];
    ha1 = [];
    hf2 = [];
    ha2 = [];
    % collect points capable of mss
    if isfield(allnss,sprintf('iid%d',iid))
        msspts = find(allnss.(['iid' num2str(iid)]));
        sslps = allnss.(['iid' num2str(iid)])(msspts);
        ss = unique(sslps);
        nss = length(ss);
        allmsspts = [];
        for iss = 1:nss
            allmsspts = union(allmsspts,msspts(sslps==ss(iss)));
            ivalpts = zeros(2*nvar,npts);
            xeqpts = zeros(2*nvar,npts);
            eqid = zeros(2,npts);
            ivalid = zeros(2,npts);
            % perturbation for all points
            for ipt = 1:npts
                pvec = alliidpvec(ipt,:,iid);                            
                % if point not capable of mss
                if ~ismember(ipt,allmsspts)                   
                    % perturbations from ss 
                    [ivalpts,ivalid,xeqpts,eqid,hf1,ha1] = ParameterPerturbations(model,pvec,...
                        xss,ivalpts,ivalid,xeqpts,eqid,ipt,tspanf,colorSpec,opts,hf1,ha1);
                    
                    % do continuation anyways to confirm
                    clear pvec
                    pvec = alliidpvec(ipt,:,iid);
                    model.PM(ac-length(orig_saddle)) = pvec(ap);
                    if eqid(1,ipt)~=eqid(2,ipt)                        
                        fprintf('Point %d is bistable. Performing continuation on\n',ipt);
                        ieq = 0;
                        while ieq<2                      
                            [data,y,p] =...
                            execMATCONT(xeqpts(nvar*ieq+1:nvar*(ieq+1),ipt),pvec,ap,fluxg,model);
                            if ~isempty(data) && size(data.s1,1)>2
                                hbif = bifurcationPlot(data.x1,data.s1,data.f1,[4,1]); 
                                fprintf('Solution %d...',ieq+1);
                                fprintf('Figure %d\n',hbif);
                            end                            
                            ieq = ieq+1;
                        end
                    else
                        fprintf('Point %d is not bistable.\n',ipt);
                        [data,y,p] =...
                        execMATCONT(xeqpts(1:nvar,ipt),pvec,ap,fluxg,model);
                        if ~isempty(data) && size(data.s1,1)>2
                            bifurcationPlot(data.x1,data.s1,data.f1,[4,1]);  
                        end
                    end                    
                else
                    s1 =...
                    siid.(['iid' num2str(iid)]).(['pt' num2str(ipt)]).s1;
                    x1 =...
                    siid.(['iid' num2str(iid)]).(['pt' num2str(ipt)]).x1;
                    f1 =...
                    siid.(['iid' num2str(iid)]).(['pt' num2str(ipt)]).f1;
                    index =...
                    cat(1,siid.(['iid' num2str(iid)]).(['pt' num2str(ipt)]).s1.index);
                    % original bifurcation plot
%                     bifurcationPlot(x1,s1,f1,[4,2]);
                    bifurcationPlot(x1,s1,f1,[4,1],ap);
%                     bifurcationPlot(x1,s1,f1,[4,3]); 
                    
                    % perturbations from ss 
                    [ivalpts,ivalid,xeqpts,eqid,hf1,ha1] = ParameterPerturbations(model,pvec,...
                        xss,ivalpts,ivalid,xeqpts,eqid,ipt,tspanf,colorSpec,opts,hf1,ha1);
                end                
            end          
        end
        % plot points from xeqpts and ivalpts using ivalid and eqid after
        % normalization
        normeqpts = xeqpts; % ./repmat(max(xeqpts,[],2),1,size(xeqpts,2));
        normival = ivalpts; % ./repmat(max(xeqpts,[],2),1,size(ivalpts,2));
        hf1 = [];
        ha1 = [];
        hf2 = [];
        ha2 = [];
        hf3 = [];
        ha3 = [];
        plotpoints = zeros(1,size(xeqpts,2));
        Point.Marker = '.';
        Point.MarkerSize = 25;
        Point2.Marker = '.';
        Point2.MarkerSize = 25;
        for ipt = 1:size(xeqpts,2) 
            addanot.text = ['P' num2str(ipt)];
            if eqid(1,ipt)==eqid(2,ipt)
                if eqid(1,ipt)==1
                    Point.MarkerFaceColor = colorSpec{1};   
                    Point.MarkerEdgeColor = colorSpec{1}; 
                elseif eqid(1,ipt)==2
                    Point.MarkerFaceColor = colorSpec{2};   
                    Point.MarkerEdgeColor = colorSpec{2}; 
                end 
                ival1 = normival(1:nvar,ipt);
                ival2 = normeqpts(1:nvar,ipt);
                plotpoints(ipt) = 1;
            else
                if eqid(1,ipt)==1
                    Point.MarkerFaceColor = colorSpec{1};   
                    Point.MarkerEdgeColor = colorSpec{1};
                    ival1 = normival(1:nvar,ipt);
                    ival2 = normeqpts(1:nvar,ipt);
                elseif eqid(1,ipt)==2
                    Point.MarkerFaceColor = colorSpec{2};   
                    Point.MarkerEdgeColor = colorSpec{2};
                    ival1 = normival(1:nvar,ipt);
                    ival2 = normeqpts(1:nvar,ipt);
                end
                if eqid(2,ipt)==1
                    Point2.MarkerFaceColor = colorSpec{1};   
                    Point2.MarkerEdgeColor = colorSpec{1}; 
                    ival3 = normival(nvar+1:end,ipt);
                    ival4 = normeqpts(nvar+1:end,ipt);
                elseif eqid(2,ipt)==2
                    Point2.MarkerFaceColor = colorSpec{2};   
                    Point2.MarkerEdgeColor = colorSpec{2};
                    ival3 = normival(nvar+1:end,ipt);
                    ival4 = normeqpts(nvar+1:end,ipt);
                end
            end
            [hf1,ha1] =...
            FIGmssEqIvalPerturbations(ival1,ival2,2,[1 2],hf1,ha1,Point,addanot); 
            [hf2,ha2] =...
            FIGmssEqIvalPerturbations(ival1,ival2,2,[2 3],hf2,ha2,Point,addanot); 
            [hf3,ha3] =...
            FIGmssEqIvalPerturbations(ival1,ival2,2,[1 3],hf3,ha3,Point,addanot); 
            if ~plotpoints(ipt)
                [hf1,ha1] =...
                FIGmssEqIvalPerturbations(ival3,ival4,2,[1 2],hf1,ha1,Point2,addanot); 
                [hf2,ha2] =...
                FIGmssEqIvalPerturbations(ival3,ival4,2,[2 3],hf2,ha2,Point2,addanot);
                [hf3,ha3] =...
                FIGmssEqIvalPerturbations(ival3,ival4,2,[1 3],hf3,ha3,Point2,addanot);
            end                 
        end
    end
end

% print analysis of results - partially complete
% same final state for both starting states
samestates = find(eqid(1,:)==eqid(2,:));

% same as strating state
highstate = find(eqid(1,:)==ivalid(1,1));
lowstate = find(eqid(2,:)==ivalid(2,1));

% systems restricted to the high state for all ivals
samehighstate = highstate(ismember(highstate,samestates));
% systems restricted to the low state for all ivals
samelowstate = lowstate(ismember(lowstate,samestates));
% systems not restricted to either of the 2 states - resting state depends
% on ival
diffstate = setdiff(1:size(eqid,2),union(samehighstate,samelowstate));

% check ival for diffstate
% systems were a low state start gets a low state and high state start
% gets a high state - i.e. maintain bistability
bistable = all(eqid(:,diffstate)==ivalid(:,diffstate));
if any(bistable)
    % no movement/change in separatrix?
    % no change in ability to remove bistability
    % system too close to old state?   
    fprintf('Following perturbations still result in bistability: %s\n'...
            ,num2str(diffstate(bistable))); 
    fprintf('Bistable parameter sets:\n');
    bistates = diffstate(bistable);
    fprintf('kEcat \t   vFbpmax \t   vEXmax\n');
    for is = 1:length(bistates)
        fprintf('%s\n',num2str(alliidpvec(bistates(is),idp),'%4.2e\t'));
    end
else
    fprintf('No perturbation results in a bistable system\n');
end
% systems restricted to high state
fprintf('Systems restricted to the high state:\n')
fprintf('%s\n',num2str(samehighstate,'%d\t'));
fprintf('High state parameters\n');
for ih = 1:length(samehighstate)
    fprintf('%s\n',num2str(alliidpvec(samehighstate(ih),idp),'%4.2e\t'));
end
% systems restricted to low state
fprintf('Systems restricted to the low state:\n')
fprintf('%s\n',num2str(samelowstate,'%d\t'));
for il = 1:length(samelowstate)
    fprintf('%s\n',num2str(alliidpvec(samelowstate(il),idp),'%4.2e\t'));
end

%% Continuation on enzyme parameters - SI Figures or Figure 4?
% enzymecont
% enzyme parameter continuation on perturbed systems
% data for Figure 3
load('C:\Users\shyam\Documents\Courses\CHE1125Project\Results\KotteModel\VmaxVariation_IvalPerturbation_Aug30.mat');
% Finer detail - Figure 3
% load('C:\Users\shyam\Documents\Courses\CHE1125Project\Results\KotteModel\VmaxVariationAll_Aug02.mat');
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
colorSpec = chooseColors(5,{'Green','Purple','Red','Navy','HotPink'});
ac = find(strcmpi(model.mets,'ac[e]'));
npts = 1;

allxeq = zeros(length(M),npts);
allxdyn = zeros(length(M),length(tspan),npts);
allxf = zeros(length(M),npts);
allfeq = zeros(length(fluxg),npts);
allfdyn = zeros(length(fluxg),length(tspan),npts);
allpvec = pvec;

% continuation on acetate
ap = 9;
solveEquilibriumODE 
ap1 = ap;
clear ap

% conitnuation on enzyme parameters 12, 13 and 14 for default values of [1
% 1 1]'
model.PM(ac-length(orig_saddle)) = orig_saddlepar;
for ip = 1:length(idp)
    ap = idp(ip);
    solveEquilibriumODE 
end
% continuation on enzyme parameters 12, 13 and 14 for perturbed values in
% cmb
% change ipt from 1 through 4 to cycle through different types of perturbations
% see Table 1 in Methods for the types fo perturbations
hf1 = [];
hf2 = [];
hf3 = [];
ipt = 4;
while ipt<size(cmb,1)
% for ipt = 1:size(cmb,1)
    xeq = alliidxeq(:,ipt);
    pvec = alliidpvec(ipt,:);
    addanot.text = ['E' num2str(ipt)];
    for ip = 1:length(idp)
        ap = idp(ip);   
        if ap~=14
            pvec(ap) = 0.01;
        else
            pvec(ap) = 5;
        end
        % run MATCONT
        [data,y,p] = execMATCONT(xeq,pvec,ap,fluxg,model);
        if ~isempty(data) && size(data.s1,1)>2
            if ap == 12
                hf1 = bifurcationPlot(data.x1,data.s1,data.f1,[4,1],ap,hf1,addanot);
            elseif ap == 13
                hf2 = bifurcationPlot(data.x1,data.s1,data.f1,[4,1],ap,hf2,addanot);
            elseif ap == 14
                hf3 = bifurcationPlot(data.x1,data.s1,data.f1,[4,1],ap,hf3,addanot);
            end            
        end
        % save MATCONT results
        s.(['pt' num2str(ipt)]) = data;
    end
    ipt = ipt+4;
end

