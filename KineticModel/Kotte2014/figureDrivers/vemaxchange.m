% vemax bifurcation diagrams and perturbation results for design
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

% perturb saddle to get steady states
eps = 1e-4;
tspanf = 0:0.1:2000;
pival = orig_saddle+eps*[1;1;1];
[~,xeq1,~,feq1] = solveODEonly(1,pival,model,pvec,opts,tspanf);
nival = orig_saddle-eps*[1;1;1];
[~,xeq2,~,feq2] = solveODEonly(1,nival,model,pvec,opts,tspanf);
xss = [xeq1 xeq2];

idp = 6;
type = 'together';
cmb = linspace(75,100,50)';
npts = size(cmb,1);

% set acetate conentration to saddle
pvec = [K1ac,K3fdp,L3,K3pep,...
        K2pep,vemax,KeFDP,ne,acetate,d,...
        k4cat,k1cat,v3max,v2max];
    
% set parameters from cmb at idp position(s)
allpvec = repmat(pvec,npts,1);
allpvec(:,idp) = cmb;

% find equilibrium point for lowest acetate 
allpvec(:,ap) = 0.001;
model.PM(ac-length(orig_saddle)) = 0.001;    
[~,allxeqlac] = solveODEonly(npts,M,model,allpvec,opts,tspan);

% and continue on acetate
[s,mssid,nss] = setupMATCONT(allxeqlac,allpvec,ap,model,fluxg,npts,10000);

%% get boundaries of acetate bistability
vemax = allpvec(mssid,idp);
acbounds = zeros(2,length(mssid)); % [min;max];
% xbounds = zeros(nvar,2*length(mssid));
mssipt = 1;
for ipt = 1:npts
    if ismember(ipt,mssid)
        index = cat(1,s.(['pt' num2str(ipt)]).s1.index);
        x1 = s.(['pt' num2str(ipt)]).x1;
%         xcont = x1(1:nvar,index);
        pcont = x1(nvar+1:end,index);
%         xbounds(:,2*mssipt-1:2*mssipt) = xcont(:,2:end-1);
        acbounds(1,mssipt) = min(pcont(:,2:end-1));
        acbounds(2,mssipt) = max(pcont(:,2:end-1));
        mssipt = mssipt+1;
    end
end
figure
plot(acbounds(1,:),vemax,'Color','b','LineWidth',2);
hold on
plot(acbounds(2,:),vemax,'Color','b','LineWidth',2);
xlabel('acetate a.u.');
ylabel('vemax s-1');

fprintf('Summary of vemax Perturbation\n')
fprintf('\t\t k4cat \t acetate\t\n');
fprintf('Minimum Bistable \t %4.3f \t %4.3f\n',min(vemax),min(acbounds(1,:)));
fprintf('Maximum Bistable \t %4.3f \t %4.3f\n',max(vemax),max(acbounds(2,:)));

%% plot bifurcation (w/ 2 LPs) on acetate
hf1 = [];
hf2 = [];
for ipt = 1:npts
    addanot.text = ['R' num2str(ipt)];
    if ismember(ipt,mssid)
        s1 = s.(['pt' num2str(ipt)]).s1;
        x1 = s.(['pt' num2str(ipt)]).x1;
        f1 = s.(['pt' num2str(ipt)]).f1;
        flux1 = s.(['pt' num2str(ipt)]).flux;
        index = cat(1,s.(['pt' num2str(ipt)]).s1.index);
        hf1 = bifurcationPlot(x1,s1,f1,[4,1],ap,hf1,addanot);  
    end
end


    