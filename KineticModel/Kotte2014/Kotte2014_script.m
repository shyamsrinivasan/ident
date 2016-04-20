addpath(genpath('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel'));
rxfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\Kotte2014\Kotte2014.txt';
cnfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\Kotte2014\Kotte2014C.txt';

% create model structure
[FBAmodel,parameter,variable,nrxn,nmetab] = modelgen(rxfname);

% obtain conentrations from file
[mc,FBAmodel,met] = readCNCfromFile(cnfname,FBAmodel);

% remove metabolites held constant from consideration in the model
% integration phase
[model,pvec,newmc,cnstmet] =...
remove_eMets(FBAmodel,parameter,mc,[FBAmodel.Vind FBAmodel.Vex],...
{'enz1[c]','ac[e]','fdp[e]'});

% only initialize for varmets   
nvar = length(model.mets)-length(find(cnstmet));
M = newmc(1:nvar);
PM = newmc(nvar+1:end);
model.PM = PM;

pvec.d = 0.25;
allhandles = feval(@Kotte2014Ckinetics);
rhsfunc = allhandles{2};
func = @(t,x)rhsfunc(t,x,model,pvec);
dMdt = func(0,M);
allmc = [M;PM];
flux = Kotte2014Ckinetics_flux(allmc,model,pvec); 

% M(1)  = 1;      % E
% M(2)  = 0.001;   % PEP
% M(3)  = 10;   % FBP
% %parameters
% kEcat = 1;
% KEacetate = 0.1;    % or 0.02
% KFbpFBP = 0.1;
% vFbpmax = 1;
% Lfbp = 4e6;
% KFbpPEP = 0.1;
% vEXmax = 1;
% KEXPEP = 0.3;
% vemax = 1.1;        % for bifurcation analysis: 0.7:0.1:1.3
% KeFBP = 0.45;       % or 0.45
% ne = 2;             % or 2
% acetate = 0.1;      % a.u acetate
% d = 0.25;           % or 0.25 or 0.35
% pvec = [kEcat,KEacetate,...
%         KFbpFBP,vFbpmax,Lfbp,KFbpPEP,...
%         vEXmax,KEXPEP,...
%         vemax,KeFBP,ne,acetate,d];
% 
% allhandles = feval(@Kotte2014glycolysis);
% rhsfunc = allhandles{2};
% func = @(t,x)rhsfunc(t,x,kEcat,KEacetate,...
%         KFbpFBP,vFbpmax,Lfbp,KFbpPEP,...
%         vEXmax,KEXPEP,...
%         vemax,KeFBP,ne,acetate,d);
% dMdt1 = func(0,M);  
% flux1 = Kotte_glycolysisflux(M,pvec);
opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
[tout,yout] = ode45(func,0:0.1:100000,M,opts);
% calculate fluxes
fout = zeros(length(tout),4);
for it = 1:length(tout)
    fout(it,:) = Kotte2014Ckinetics_flux([yout(it,:)';PM],model,pvec);
end
plot(tout,yout);

plotKotteVariables(tout,yout,1);