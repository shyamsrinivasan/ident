% find points on the stability boundary for Kotte model
% Chiang, Hirsch and Wu, 1988 (refered from Khalil, Nonlinear Systems)
% step 1. find all equilibrium points
% build stoichioemtrc matrices
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
k1cat = 1;
K1ac = 0.1;    % or 0.02
K3fdp = 0.1;
v3max = 1;
L3 = 4e6;
K3pep = 0.1;
v2max = 1;
K2pep = 0.3;
vemax = 1.1;        % for bifurcation analysis: 0.7:0.1:1.3
KeFDP = 0.45;       % or 0.45
ne = 2;             % or 2
acetate = 0.1;      % a.u acetate
d = 0.25;           % or 0.25 or 0.35
k4cat = 0.2;
pvec = [K1ac,K3fdp,L3,K3pep,...
        K2pep,vemax,KeFDP,ne,acetate,d,...
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
[orig_saddle,orig_saddlepar,saddleid] = getsaddlenode(data.s1,data.x1,5e-3);
pvec(ap) = orig_saddlepar;
model.PM(ac-length(orig_saddle)) = orig_saddlepar;
saddleflux = Kotte_givenFlux([orig_saddle;model.PM],pvec,model);
% get eigen values and eigen vectors of the saddle node
[~,alleig,alleigw] = KotteStabilityInfo(orig_saddle',model,pvec);

% perturb saddle to get steady states
eps = 1e-4;
tspanf = 0:0.1:2000;
pival = orig_saddle+eps*[1;1;1];
[~,xeq1,~,feq1] = solveODEonly(1,pival,model,pvec,opts,tspanf);
nival = orig_saddle-eps*[1;1;1];
[~,xeq2,~,feq2] = solveODEonly(1,nival,model,pvec,opts,tspanf);
xss = [xeq1 xeq2];

%% pick points on the boundary of region of attraction - unstable eigen vectors
figure
hold on
itermax = 100;
eps = 1e-3;
alpha = 0.8;
[ip_interUS,status] =...
getUSintsecpt(orig_saddle,alleigw(:,3),model,pvec,opts,-1,eps,alpha,itermax);
for ip = 1:size(ip_interUS,2)
    xfwd = solveODEonly(1,ip_interUS(:,ip),model,pvec,opts,[0,2000]);
    xrev = solveODEonly(1,ip_interUS(:,ip),model,pvec,opts,[0,-40]);
%     plot(ip_interUS(1,ip),ip_interUS(2,ip),'Marker','.','MarkerSize',10,'Color','b');
%     plot(xfwd(1,:),xfwd(2,:),'k');
%     plot(xrev(1,:),xrev(2,:),'r');
    plot3(ip_interUS(1,ip),ip_interUS(2,ip),ip_interUS(3,ip),'Marker','.','MarkerSize',10,'Color','b');
    plot3(xfwd(1,:),xfwd(2,:),xfwd(3,:),'k');
    drawnow
    plot3(xrev(1,:),xrev(2,:),xrev(3,:),'r');
    drawnow
end

%% pick points on the boundary of region of attraction - stable eigen vectors
% figure
% hold on
itermax = 5;
eps = 1e-3;
alpha = 0.8;
for iw = 1:2
    [ip_interS,status] =...
    getSintsecpt(orig_saddle,alleigw(:,iw),model,pvec,opts,[1 2],eps,alpha,itermax);
    for ip = 1:size(ip_interUS,2)
        xfwd = solveODEonly(1,ip_interS(:,ip),model,pvec,opts,[0,2000]);
        xrev = solveODEonly(1,ip_interS(:,ip),model,pvec,opts,[0,-40]);
    %     plot(ip_interS(1,ip),ip_interS(2,ip),'Marker','.','MarkerSize',10,'Color','b');
    %     plot(xfwd(1,:),xfwd(2,:),'k');
    %     plot(xrev(1,:),xrev(2,:),'r');
        plot3(ip_interS(1,ip),ip_interS(2,ip),ip_interS(3,ip),'Marker','.','MarkerSize',10,'Color','b');
        plot3(xfwd(1,:),xfwd(2,:),xfwd(3,:),'k');
        drawnow
        plot3(xrev(1,:),xrev(2,:),xrev(3,:),'r');
        drawnow
    end
end