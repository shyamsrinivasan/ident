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

idp = 4;
type = 'together';
cmb = linspace(0,0.5,30)';
npts = size(cmb,1);

% set acetate conentration to saddle
pvec = [K1ac,K3fdp,L3,K3pep,...
        K2pep,vemax,KeFDP,ne,acetate,d,...
        k4cat,k1cat,v3max,v2max];
pvec(3) = 3e6;
% pvec(ap) = orig_saddlepar;
% model.PM(ac-length(orig_saddle)) = orig_saddlepar;

% set parameters from cmb at idp position(s)
allpvec = repmat(pvec,npts,1);
allpvec(:,idp) = cmb;

% find equilibrium point for lowest acetate 
allpvec(:,ap) = 0.01;
model.PM(ac-length(orig_saddle)) = 0.01;    
[~,allxeqlac] = solveODEonly(npts,M,model,allpvec,opts,tspan);

% and continue on acetate
[s,mssid,nss] = setupMATCONT(allxeqlac,allpvec,ap,model,fluxg,npts,1500);

%% % plot all available steady state pep concentrations for a 
% given parameter V4max against all available fluxes
% choose ipt = 1, 5, 11, 15 and 17 (all points 1:17 are bistable)
figure
hold on
for ipt = 3:7 % npts
    if ismember(ipt,mssid)
        bifpts = cat(1,s.(['pt' num2str(ipt)]).s1.index);
        for i = 1:3
            x1 = s.(['pt' num2str(ipt)]).x1(:,bifpts(i)+1:bifpts(i+1)-1);
            flux = s.(['pt' num2str(ipt)]).flux(:,bifpts(i)+1:bifpts(i+1)-1);
            eig = s.(['pt' num2str(ipt)]).f1(:,bifpts(i)+1:bifpts(i+1)-1);            
            if any(real(eig)>0)
                style = '--';
            else
                style = '-';
            end
            hl1 = line(x1(1,:),flux(1,:),'Color',colorSpec{1},'LineStyle',style);
            hl2 = line(x1(1,:),flux(4,:),'Color',colorSpec{2},'LineStyle',style);
            hl3 = line(x1(1,:),flux(5,:),'Color',colorSpec{3},'LineStyle',style);
        end
        xlps = s.(['pt' num2str(ipt)]).x1(:,bifpts);
        flps = s.(['pt' num2str(ipt)]).flux(:,bifpts);
        line(xlps(1,[2,3]),...
             flps(1,[2,3]),'Color','r','Marker','.','LineStyle','none');
        line(xlps(1,[2,3]),...
             flps(4,[2,3]),'Color','r','Marker','.','LineStyle','none');
        line(xlps(1,[2,3]),...
             flps(5,[2,3]),'Color','r','Marker','.','LineStyle','none');
        legend([hl1 hl3 hl2],'v1','v3','v2');     
        drawnow
    end
end

%% % plot all available steady state pep concentrations for a 
% given parameter V4max against all available concentrations
% choose ipt = 1, 5, 11, 15 and 17
figure
hold on
for ipt = 3:7 % npts
    if ismember(ipt,mssid)
        bifpts = cat(1,s.(['pt' num2str(ipt)]).s1.index);
        for i = 1:3
            x1 = s.(['pt' num2str(ipt)]).x1(:,bifpts(i)+1:bifpts(i+1)-1);
            flux = s.(['pt' num2str(ipt)]).flux(:,bifpts(i)+1:bifpts(i+1)-1);
            eig = s.(['pt' num2str(ipt)]).f1(:,bifpts(i)+1:bifpts(i+1)-1);            
            if any(real(eig)>0)
                style = '--';
            else
                style = '-';
            end
            hl1 = line(x1(1,:),x1(2,:),'Color',colorSpec{1},'LineStyle',style);
            hl2 = line(x1(1,:),x1(3,:),'Color',colorSpec{2},'LineStyle',style);           
        end
        xlps = s.(['pt' num2str(ipt)]).x1(:,bifpts);
        flps = s.(['pt' num2str(ipt)]).flux(:,bifpts);
        line(xlps(1,[2,3]),...
             xlps(2,[2,3]),'Color','k','Marker','.','LineStyle','none');
        line(xlps(1,[2,3]),...
             xlps(3,[2,3]),'Color','k','Marker','.','LineStyle','none');        
        legend([hl1 hl2],'fdp','E');     
        drawnow
    end
end
