% get different stable maniforlds from Kotte model using information from
% Lyons et al., 2014, Int. J. Bifurcation Chaos

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

% all manifold in 3D
% tspanr = [0,-30;0 -12];
% NumericalSeparatrix(model,pvec,opts,ap,data.s1,data.x1,[xeq1 xeq2],'all',tspanr,3,5e-3);

% Lyons et al., 2014 code
% mypath = 'C:\Users\shyam\Documents\MATLAB\CCFM_manifolds';
% addpath(strcat(mypath,'\CCFM_manifolds\functions\'));

% set viewpoint for 3D plots
azimuth = 52; elevation = 10;

% color_hash for plotting
color_hash = {'k.','go','ko'};

% compute invariant manifold via numerical integration
ti = 0.0; h = -0.5; tf = -25;%-25;%integration time data
time = ti:h:tf; %set time array
save_time_start_index = 3;

eqpts = [saddle'; % 2D stable 
       xeq1';   % spiral sink (3D stable)
       xeq2'];  % sink (3D stable)

npoints = 500;
range = [saddlepar 2];
contdir = 1;
options = optimoptions('fsolve','Display','off','TolFun',1e-10,'TolX',1e-10);
fixed_points = kotte_branches(npoints,range,contdir,eqpts,model,pvec,ap,options);   
[type,alleig] = KotteStabilityInfo(eqpts,model,pvec);      

% compute Jacobian, eigenvalues, eigenvector at saddle points with 2D stable
% manifolds
[~,eig,eigvec] = getKotteJacobian(eqpts(1,:)',pvec,model);

% obtain the 2 eigenvectors of 2D linear eigenspace
stableeigvec = eigvec(:,eig<0);

% circle parameters
points = 401;
radius = 0.01;

% 2D stable manifold 
int_time_1D_manifolds = 0:-0.05:-15;
time_1D = length(int_time_1D_manifolds);
one_dim_manifolds_array_x = zeros(2*time_1D,6);
one_dim_manifolds_array_y = zeros(2*time_1D,6);
one_dim_manifolds_array_z = zeros(2*time_1D,6);

% one single saddle node with 2D stable manifold
delta_x = 0.01;
step = 0.005;
pvec(ap) = saddlepar;
model.PM(ac-length(saddle)) = saddlepar;

% input functions
givenMfsolve = @(x)Kotte_givenNLAE(x,model,pvec);
ode_system_fsolve = 'givenMfsolve';
givenModel = @(t,x)KotteODE(t,x,model,pvec);
ode_system = 'givenModel';

% compute eigenvalues for 2D stable manifolds with integrating backwards
% in time beginning with a circle of initial conditions with radius "r"
% and "points" points in the 2D space spanned by the two eigenvectors
% at the equilibrium "fixed_points" with exactly 2 eigenvalues with
% negative real part

% % obtain coordinates of circle with radius r in (x1,x2) plane
[x1,x2] = getplanarcircle(points,radius);
% 
% % perform linear mapping of unit circle onto plane in R3 spanned by W1,W2
xnew = manifoldlinearmapping(x1,x2,stableeigvec(:,1),stableeigvec(:,2));
% 
% % translate mapping to saddle point
xnew = xnew + repmat(saddle,1,size(xnew,2));
xnew = xnew';
% 
tspanr = 0:-.05:-25;
[x,y,z] = get2Dmanifoldpoints(xnew,model,pvec,tspanr,opts);

% [xnew,ynew,znew] = removeredundantpoints(x,y,z,0.01);

[xnew,ynew,znew] = redundant_point_filter(x,y,z,0.01);

figure(4); hold on;
Delaunay_special_plot(real(xnew),real(ynew),real(znew),0.05);


% trim dynr
% txdynr = allxdynr;
% txdynr(:,txdynr(1,:)>100) = [];
% 

% 
% 
% % filter out extra points within small tolerance of each other
% norm_tol = 0.1;%0.01;
% [x,y,z] = redundant_point_filter(allxdynr(1,:),allxdynr(2,:),allxdynr(3,:),norm_tol);
% 
% figure(3); hold on;
% Delaunay_special_plot(x,y,z,0.05);






                       
                       