% generate equilibrium solution and model for Kotte model
if ~exist('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\Kotte2014\Kotte2014.txt')
    status = 2;
    fprintf('\nLinux System\n');
else 
    status = 1;
    fprintf('\nWindows System\n');
end
if status == 1
    rxfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\Kotte2014\Kotte2014.txt';
    cnfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\Kotte2014\Kotte2014C.txt';
elseif status == 2
    rxfname = '/home/shyam/Documents/MATLAB/Code/KineticModel/Kotte2014/Kotte2014.txt';
    cnfname = '/home/shyam/Documents/MATLAB/Code/KineticModel/Kotte2014/Kotte2014C.txt';
end

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
givenMfsolve = @(x)Kotte_givenNLAE(x,model,pvec);
fluxg = Kotte_givenFlux([M;model.PM],pvec,model);
dMdtg = givenModel(0,M);

tspan = 0:0.1:500;    
npts = 1;
allpvec = pvec;  

% run equilibrium solution followed by MATCONT
opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
ac = find(strcmpi(model.mets,'ac[e]'));
allxeq = zeros(length(M),npts);
allxdyn = zeros(length(M),length(tspan),npts);
allxf = zeros(length(M),npts);
allfeq = zeros(length(fluxg),npts);
allfdyn = zeros(length(fluxg),length(tspan),npts);
ap = 9;
solveEquilibriumODE    

% get saddle node
[saddle,saddlepar] = getsaddlenode(data.s1,data.x1,5e-3);
pvec(ap) = saddlepar;
model.PM(ac-length(saddle)) = saddlepar;

% calculate eig val for all points on continuation curve data.x1
alleig = zeros(3,size(data.x1,2));
for ipts = 1:size(data.x1,2)
    pvec(ap) = data.x1(end,ipts);
    model.PM(ac-length(saddle)) = data.x1(end,ipts);
    [~,alleig(:,ipts)] = KotteStabilityInfo(data.x1(1:3,ipts)',model,pvec);      
end

% perturb saddle to get steady states
pvec(ap) = saddlepar;
model.PM(ac-length(saddle)) = saddlepar;
eps = 1e-4;
tspanf = 0:0.1:2000;
pival = saddle+eps*[1;1;1];
[~,xeq1] = solveODEonly(1,pival,model,pvec,opts,tspanf);

nival = saddle-eps*[1;1;1];
[~,xeq2] = solveODEonly(1,nival,model,pvec,opts,tspanf);

[~,eigval,w] = getKotteJacobian(saddle,pvec,model);

% all manifolds in 2D
% tspanr = [0,-30;0 -12];
% NumericalSeparatrix(model,pvec,opts,ap,data.s1,data.x1,[xeq1 xeq2],'all',tspanr,2,5e-3);

% unstable manifold in 3D
tspanr = [0:-.1:-30]; % 0 -12];
[hfig,xmafold] =...
NumericalSeparatrix(model,pvec,opts,ap,data.s1,data.x1,[xeq1 xeq2],...
                    'unstable',tspanr,3,5e-3);
Perturb2DstableManifold(xmafold(1:3,:)',saddle,model,pvec,opts,tspanf,...
                        hfig,[1 2 3],10); 