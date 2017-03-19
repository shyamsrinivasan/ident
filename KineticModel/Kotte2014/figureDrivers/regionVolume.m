% get volume of region of attraction based on MC simulation of initial
% values and checking final equilibrium points
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

% get unstable manifold to determine the divide between regions of
% attraction
tspanr = 0:-.1:-30;
tspanf = 0:0.1:2000;
eps = 1e-4;
[~,eig,w] = getKotteJacobian(orig_saddle,pvec,model);
[xWus,xeq] = calc1DWus(orig_saddle,w,eig,model,pvec,opts,tspanr,tspanf,eps);
% chop xWus
[~,nzid,~] = find(xWus~=0,1,'last');
relWus = real(xWus(:,1:nzid));
[x,y,z] = chopvals(relWus(1,:),relWus(2,:),relWus(3,:),[10 10 10]);

%% get region of attraction by perturbing around the steady state and 
% expanding region of perturbation
% radius where no initial point produces the 2nd steady state = rAss1 =
% [1.22 1.22 1.22]'
rndivals = get3Dsphere(xeq1',1.22,50000);
save('regionVoluemSamplevals');
% integrate
options = [];
% hfig = figure;
% [ppival,npival,allxeq,ssid,allfeq] =...
% getPequilibrium(rndivals,model,pvec,options,opts,tspanf);
% plotrndivals(rndivals,ssid,allxeq,[1 2 3],2,hfig,[]);
% plotrndivals(rndivals,ssid,allxeq,[1 2],2,hfig,[]);
% plotrndivals(rndivals,ssid,allxeq,[2 3],2,hfig,[]);
% plotrndivals(rndivals,ssid,allxeq,[2 3],2,hfig,[]);


%% get random initial values for all 3 variables and calculate steady state
% rndivals = randomivals([0 5;0 5;0 5],10000);
% % integrate
% options = [];
% hfig = figure;
% [ppival,npival,allxeq,ssid,allfeq] =...
% getPequilibrium(rndivals,model,pvec,options,opts,tspanf);
% plotrndivals(rndivals,ssid,allxeq,[1 2],2,hfig,ha)
%% get nipts random initial values and corresponding equilibrium points
% through perturbation of relWus (x,y,z)
% relWus = [x;y;z];
% newivals = perturbTrajectory(relWus);

%% integrate from newivals
% options = optimoptions('fsolve','Display','final-detailed',...
%                        'TolFun',1e-16,...
%                        'TolX',1e-12,...
%                        'MaxFunEvals',1000000,...
%                        'MaxIter',50000);
% hfig = figure;
% [allxeq,ssid,allfeq] =...
% getPequilibrium(newivals,model,pvec,options,opts,tspanf,hfig,[1 2 3]);
