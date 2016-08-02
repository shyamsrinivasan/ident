% enzyme parameter continuation

% data from perturbation continuation on acetate
% Figure 3
load('C:\Users\shyam\Documents\Courses\CHE1125Project\Results\KotteModel\VmaxVariationAllPerturbations_July29.mat');

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
% continuation on acetate
ap = 9;
solveEquilibriumODE 
% conitnuation on enzyme parameters 12, 13 and 14 for default values of [1
% 1 1]'
for ip = 1:length(idp)
    ap = idp(ip);
    solveEquilibriumODE 
end
% continuation on enzyme parameters 12, 13 and 14 for perturbed values in
% cmb
hf1 = [];
hf2 = [];
hf3 = [];
for ipt = 1:size(cmb,1)
    xeq = alliidxeq(:,ipt);
    pvec = alliidpvec(ipt,:);
    addanot.text = ['P' num2str(ipt)];
    for ip = 1:length(idp)
        ap = idp(ip);        
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
end
    
