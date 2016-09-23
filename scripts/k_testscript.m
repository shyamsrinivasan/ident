clc
% clear all
addpath(genpath('C:\Users\shyam\Documents\Courses\CHE 1125 Project\Kinetic Model'));
% load('C:\Users\shyam\Documents\Courses\CHE 1125 Project\Kinetic Model\kmodel_pname.mat');
rxfname = 'C:\Users\shyam\Documents\Courses\CHE 1125 Project\Kinetic Model\rxnlist.txt';
[FBAmodel,parameter,variable,nrxn,nmetab] = modelgen(rxfname);
% % % sample metabolites 
% variable.MC = sampleMetabolites(FBAmodel);
variable.MC = sampleM_1(FBAmodel);
nmodels = 3;
[ensb] = build_ensemble(nmodels,FBAmodel,parameter,variable.MC);
%Model initialization
FBAmodel.uptake_flux = 20;%mmole/h ->still confused!
FBAmodel.gmax = 0.8;%h-1
% Check Model Stability of ensemble
% Use jacobians and their Eigen values
% ODE solver parameters
SolverOptions.scalAbsTol = 1e-8;
SolverOptions.RelTol = 1e-8;
SolverOptions.MaxIter = 1000;    
SolverOptions.MaxDataPoints = 200; 
SolverOptions.tout = 0.01;
SolverOptions.tmax = 10;%s
%define perturbations
pertb.pertb1.enzname = 'Protein2';
pertb.pertb1.change = 'increase';
pertb.pertb1.percent = 0;
pertb.pertb2.enzname = 'Protein2';
pertb.pertb2.change = 'increase';
pertb.pertb2.percent = 100;
pertb.pertb3.enzname = 'Protein2';
pertb.pertb3.change = 'decrease';
pertb.pertb3.percent = 100;
pertb.pertb4.enzname = 'Protein5';
pertb.pertb4.change = 'increase';
pertb.pertb4.percent = 100;
[inSolution,enSolution] =  solveEnsemble(ensb,FBAmodel,variable,pertb,SolverOptions);
%Compare Solution from different models


%ADMAT
Y = enSolution.model1.y(:,end);
nvar = size(Y,1);
obj = deriv(Y,eye(nvar));
pmeter = ensb.model1;
Jy = ODEmodel(0,obj,FBAmodel,pmeter);


Extra = struct();
Extra.model = model;
Extra.pmeter = ensb.model1;
[dXdt,Jx] = feval(ADfun('ODEmodel_ADMAT',nvar),Y,Extra);


%Estimate/sample unknown quantities
%Unknown Quantities: delG, Keq, Metabolite Concentrations (MC)
%             MC - sign(Vnet) = sign(RTlnKeq/pi(MC^n))

%Sampling enzyme saturation sigma
%Identify reactions with inhibitors and consider to act as competetive only
%Allosteric inhibition - Inhibitors/activators bind to different sites than
%substrates
%NADH.ATP,ADP,NAD+,NADP,NADPH - all bind as ping pong mechanisms => single
%active site


%[Jthermo, Jsat, Jreg] = calcJacobian(kinmodel);
%Flux using convinience kinetics
% MC = model.MC;
% [initflux] = ConvinienceKinetics(model,biomassind,MC);
%runOpt(kinmodel,initflux);
%v = ConvinienceKinetics(kinmodel,Conc,Parameters);
clc
% addpath(genpath('C:\Users\shyam\Documents\MATLAB\sundialsTB'));
% metabfname = 'C:\Users\shyam\Documents\Courses\CHE 1125 Project\Metabolic Model\metabolite_conc.txt';
% [kmodel,elasticity,tnspecies,tnsteps] = KModelGen(metabfname,fbamodel,nrxn,nmetab);