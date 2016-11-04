% get different stable maniforlds from Kotte model using information from
% Lyons et al., 2014, Int. J. Bifurcation Chaos

% generate equilibrium solution and model for Kotte model
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

% perturb saddle to get steady states
eps = 1e-4;
tspanf = 0:0.1:2000;
pival = saddle+eps*[1;1;1];
[~,xeq1] = solveODEonly(1,pival,model,pvec,opts,tspanf);

nival = saddle-eps*[1;1;1];
[~,xeq2] = solveODEonly(1,nival,model,pvec,opts,tspanf);

% Lyons et al., 2014 code
mypath = 'C:\Users\shyam\Documents\MATLAB\CCFM_manifolds';
addpath(strcat(mypath,'\CCFM_manifolds\functions\'));

%initialize parameters
global c lambda;
c = 3.0;

% solution tolerance (f_solve)
sol_tol = 1.0e-8;

% set viewpoint for 3D plots
azimuth = 52; elevation = 10;

% color_hash for plotting
color_hash = {'k.','go','ko'};

% compute invariant manifold via numerical integration
ti = 0.0; h = -0.25; tf = -25;%-25;%integration time data
time = ti:h:tf; %set time array
save_time_start_index = 3;

eqpts = [saddle'; % 2D stable 
       xeq1';   % spiral sink (3D stable)
       xeq2'];  % sink (3D stable)
                       
[type,alleig] = KotteStabilityInfo(eqpts,model,pvec);       

% circle parameters
points = 401;
radius = 0.05;

% 2D stable manifold 
int_time_1D_manifolds = 0:-0.005:-15;
time_1D = length(int_time_1D_manifolds);
one_dim_manifolds_array_x = zeros(2*time_1D,6);
one_dim_manifolds_array_y = zeros(2*time_1D,6);
one_dim_manifolds_array_z = zeros(2*time_1D,6);

% one single saddle node with 2D stable manifold
delta_x = 0.01;
step = 0.005;
pvec(ap) = 0.01;
model.PM(ac-length(saddle)) = 0.01;

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

% compute Jacobian, eigenvalues, eigenvector at saddle points with 2D stable
% manifolds
[~,eig,eigvec] = getKotteJacobian(eqpts(1,:)',pvec,model);

% obtain the 2 eigenvectors of 2D linear eigenspace
eigvec_stable = eigvec(:,eig<0);

% obtain coordinates of circle with radius r in (x1,x2) plane
[x1,x2] = planar_circle(points,radius);

% perform linear mapping of unit circle onto plane in R3 spanned by W1,W2
[x1_new,x2_new,x3_new] =...
curve_map_R2_R3(x1,x2,eigvec_stable(:,1),eigvec_stable(:,2),points);

%translate mapping to saddle point
x1_new = x1_new + saddle(1);
x2_new = x2_new + saddle(2);
x3_new = x3_new + saddle(3);
x_new = [x1_new' x2_new' x3_new'];

[x,y,z] =...
manifold_compute_circle(saddle',time,save_time_start_index,delta_x,step,ode_system,x_new,1);



                       
                       