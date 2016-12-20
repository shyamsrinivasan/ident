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

[xchop,ychop,zchop,r] = chopvals(x,y,z,5);

[xnew,ynew,znew] = removeredundantpoints(real(xchop),real(ychop),real(zchop),0.05);

[xxnew,yynew,zznew] = removeredundantpoints(real(xchop),real(ychop),real(zchop),0.01);

Manifold2DPlot(real(xnew),real(ynew),real(znew));


%% Separate out segments of manifold surface
id2 = ynew<2.6065;
id2 = setdiff(1:length(ynew),find(id2));

idd2 = yynew<2.6065;
idd2 = setdiff(1:length(yynew),find(idd2));

x1 = xnew;xx1 = xxnew;
y1 = ynew;yy1 = yynew;
z1 = znew;zz1 = zznew;

x1(id2) = [];y1(id2) = [];z1(id2) = [];
xx1(idd2) = [];yy1(idd2) = [];zz1(idd2) = [];

id1 = y1>1.185;
idd1 = yy1>1.1875;

% surface 1
xs1 = x1(id1);ys1 = y1(id1);zs1 = z1(id1);
xxs1 = xx1(idd1);yys1 = yy1(idd1);zzs1 = zz1(idd1);

x1(id1) = [];y1(id1) = [];z1(id1) = [];
xx1(idd1) = [];yy1(idd1) = [];zz1(idd1) = [];

idd3 = yy1>0.49155;

xxs23 = xx1(idd3);yys23 = yy1(idd3);zzs23 = zz1(idd3);

idd4 = zzs23<0.6103;
idd5 = setdiff(1:length(xxs23),find(idd4));

xxs2 = xxs23(idd5);yys2 = yys23(idd5);zzs2 = zzs23(idd5);
xxs3 = xxs23(idd4);yys3 = yys23(idd4);zzs3 = zzs23(idd4);

idd5 = yys3>1.092;
xxs4 = xxs3(idd5);yys4 = yys3(idd5);zzs4 = zzs3(idd5);

% surface 2
xxs3(idd5) = [];yys3(idd5) = [];zzs3(idd5) = [];
xxs2 = [xxs3 xxs2];yys2 = [yys3 yys2];zzs2 = [zzs3 zzs2];
clear xxs3 yys3 zzs3

% surface 3
xxs1 = [xxs1 xxs4];yys1 = [yys1 yys4];zzs1 = [zzs1 zzs4];
clear xxs4 yys4 zzs4

xx1(idd3) = [];yy1(idd3) = [];zz1(idd3) = [];

idd6 = zz1<2.000&xx1<0.9505;
xxs3 = xx1(idd6);yys3 = yy1(idd6);zzs3 = zz1(idd6);

idd7 = setdiff(1:length(xx1),find(idd6));
xxs4 = xx1(idd7);yys4 = yy1(idd7);zzs4 = zz1(idd7);

idd8 = xxs4<0.3862&zzs4<2.395;
xxs5 = xxs4(idd8);yys5 = yys4(idd8);zzs5 = zzs4(idd8);
xxs3 = [xxs3 xxs5];yys3 = [yys3 yys5];zzs3 = [zzs3 zzs5];

xxs4(idd8) = [];yys4(idd8) = [];zzs4(idd8) = [];

idd9 = xxs3<0.3273&zzs3<0.5410;
xxs6 = xxs3(idd9);yys6 = yys3(idd9);zzs6 = zzs3(idd9);

xxs3(idd9) = [];yys3(idd9) = [];zzs3(idd9) = [];

xxs7 = xxs3(zzs3>1.241);yys7 = yys3(zzs3>1.241);zzs7 = zzs3(zzs3>1.241);
xxs8 = xxs3(zzs3<1.241);yys8 = yys3(zzs3<1.241);zzs8 = zzs3(zzs3<1.241);

%% Final Figure
hfig = Manifold2DPlot(real(xxs1),real(yys1),real(zzs1),[],[1 0 0]);
Manifold2DPlot(real(xxs2),real(yys2),real(zzs2),hfig,[0 1 0]);
Manifold2DPlot(real(xxs4),real(yys4),real(zzs4),hfig,[0 0 0]);
% Manifold2DPlot(real(xxs5),real(yys5),real(zzs5),hfig,[0.5 0.5 0.5]);
Manifold2DPlot(real(xxs6),real(yys6),real(zzs6),hfig,[0 0 1]);
Manifold2DPlot(real(xxs7),real(yys7),real(zzs7),hfig,[0 0 1]);
Manifold2DPlot(real(xxs8),real(yys8),real(zzs8),hfig,[0 0 1]);
plot3(saddle(1),saddle(2),saddle(3),'Color','r','Marker','.','MarkerSize',30);






% [xnew,ynew,znew] = redundant_point_filter(x,y,z,0.01);


% figure; hold on;
% Delaunay_special_plot(real(xnew),real(ynew),real(znew),0.05);


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

% xnew coordinates for 2 surfaces
% surface 1
% i1 = [1;1294;1507;1717;1928;2137;2345;2550;2753;2954;3079;3204;3332;...
%     3458;3584;3708;3833;3958;4082;4204;4325;4446;4566;4776;4985;5193;...
%     5399;5604;5808;6010;6211;6411;6610;6808;7005;7200;7394;7587;7779;...
%     7970;8161;8350;8538;8724;8909;9094;9278;9461;9643;9824;10004;10183;10362];
% i2 = [1178;1391;1605;1815;2025;2234;2440;2644;2846;3047;3172;3296;3423;...
%     3549;3674;3798;3924;4049;4171;4293;4414;4535;4637;4847;5056;5264;...
%     5470;5675;5878;6080;6281;6481;6680;6878;7074;7269;7463;7656;7848;...
%     8039;8230;8419;8606;8792;8977;9162;9346;9529;9712;9892;10072;10251;10430];

% surface 2
% ii1
% ii2







                       
                       