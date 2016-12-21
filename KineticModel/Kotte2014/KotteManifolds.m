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

% unstable manifold in 3D
tspanr = [0,-30]; % 0 -12];
hfig =...
NumericalSeparatrix(model,pvec,opts,ap,data.s1,data.x1,[xeq1 xeq2],'unstable',tspanr,3,5e-3);

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

%%
% compute eigenvalues for 2D stable manifolds with integrating backwards
% in time beginning with a circle of initial conditions with radius "r"
% and "points" points in the 2D space spanned by the two eigenvectors
% at the equilibrium "fixed_points" with exactly 2 eigenvalues with
% negative real part
% circle parameters
points = 801;
radius = 0.01;

% % obtain coordinates of circle with radius r in (x1,x2) plane
[x1,x2] = getplanarcircle(points,radius);
% 
% % perform linear mapping of unit circle onto plane in R3 spanned by W1,W2
circlenew = manifoldlinearmapping(x1,x2,stableeigvec(:,1),stableeigvec(:,2));
% 
% % translate mapping to saddle point
circlenew = circlenew + repmat(saddle,1,size(circlenew,2));
circlenew = circlenew';
% 
tspanr = 0:-.05:-25;
[x,y,z,dynr] = get2Dmanifoldpoints(circlenew,model,pvec,tspanr,opts);
[xchop,ychop,zchop,r] = chopvals(x,y,z,[2.3 1.6 1.6]);
% [xnew,ynew,znew] = removeredundantpoints(real(xchop),real(ychop),real(zchop),0.01);
% Manifold2DPlot(real(xnew),real(ynew),real(znew));

%% Separate out segments of manifold surface
% working with xchop
% backup 
x1 = real(xchop);y1 = real(ychop);z1 = real(zchop);
id2 = y1<1.1885;
id1 = setdiff(1:length(x1),find(id2));

xs1 = x1(id1);
ys1 = y1(id1);
zs1 = z1(id1);
x1(id1) = [];y1(id1) = [];z1(id1) = [];

% x1(id1) = [];y1(id1) = [];z1(id1) = [];
% [x1,y1,z1] = removeredundantpoints(real(x1),real(y1),real(z1),0.01);

id3 = x1<1.325&y1>0.485&z1<2.6065;
xs2 = x1(id3);ys2 = y1(id3);zs2 = z1(id3);
x1(id3) = [];y1(id3) = [];z1(id3) = [];
% 
% % small segmet in xs2
id3 = xs2>.2491 & xs2<.4319 & ys2>1.0905 & ys2<1.195 & zs2>.5164 & zs2<.5796;
xs1 = [xs1 xs2(id3)];ys1 = [ys1 ys2(id3)];zs1 = [zs1 zs2(id3)];
[xs1,ys1,zs1] = removeredundantpoints(xs1,ys1,zs1,0.02);
xs2(id3) = [];ys2(id3) = [];zs2(id3) = [];

id4 = x1<0.687 & z1<2.4;
xs3 = x1(id4);ys3 = y1(id4);zs3 = z1(id4);
x1(id4) = [];y1(id4) = [];z1(id4) = [];
% divide xs3 into 2 parts
id5 = ys3>0.06834;
id6 = setdiff(1:length(ys3),find(id5));
xs31 = xs3(id5);ys31 = ys3(id5);zs31 = zs3(id5);
xs32 = xs3(id6);ys32 = ys3(id6);zs32 = zs3(id6);

hfig = Manifold2DPlot(xs1,ys1,zs1,hfig);
Manifold2DPlot(xs2,ys2,zs2,hfig);
Manifold2DPlot(xs31,ys31,zs31,hfig);
Manifold2DPlot(xs32,ys32,zs32,hfig);
Manifold2DPlot(x1,y1,z1,hfig);
% plot3(saddle(1),saddle(2),saddle(3),'Color','r','Marker','.','MarkerSize',30);

% 
% % foil on surface xs2
% id5 = xs2<1.325 & ys2>0.7312 & ys2<1.189 & zs2>0.6305 & zs2<2.606;
% xs3 = xs2(id5);ys3 = ys2(id5);zs3 = zs2(id5);
% xs2(id5) = [];ys2(id5) = [];zs2(id5) = [];

% id6 = xs2>0.0004287 & xs2<0.3842 & ys2>0.7312 & ys2<0.9922 & zs2>0.5523 & zs2<0.61375;
% xs3 = [xs3 xs2(id6)];ys3 = [ys3 ys2(id6)];zs3 = [zs3 zs2(id6)];
% xs2(id6) = [];ys2(id6) = [];zs2(id6) = [];                    
                       