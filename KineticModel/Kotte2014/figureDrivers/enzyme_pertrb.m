% enzyme perturbation script - Figure 7
% changes to enzyme expression levels to determine new steady states by
% simulating from either of the 2 initial steady states

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
        % perturbations from ss 
        [ivalpts,ivalid,xeqpts,eqid,hf1,ha1] = ParameterPerturbations(model,pvec,...
            xss,ivalpts,ivalid,xeqpts,eqid,ipt,tspanf,colorSpec,opts,hf1,ha1);                          
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

% plotting section
% run from either saved data or directly following above simulation
npts = size(alliidpvec,1);
nvar = size(alliidxeq,1);
ndp = size(alliidpvec,3);
tout = tspan;
colorSpec = chooseColors(4,{'Navy','HotPink','Red','Orange'});

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
                else                    
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
