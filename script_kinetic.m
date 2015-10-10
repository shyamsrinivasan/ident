clc

addpath(genpath('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel'));
rxfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\ecoliN1.txt';

%create model structure
[FBAmodel,parameter,variable,nrxn,nmetab] = modelgen(rxfname);

%add rxn
rxn.equation = 'atp[c] <==>';
FBAmodel = addRxn(FBAmodel,rxn);

%assign initial fluxes and calculate FBA fluxes for direction
FBAmodel = FBAfluxes(FBAmodel,'pfba');

% [x,xFlag,xmax,xmin] = getiConEstimate(FBAmodel);
% FBAmodel.x = x;

setupMetLP_2(FBAmodel)

% x = initialsample(FBAmodel);

%samwple initial metabolite concentrations for estimating kientic parameters
[mc,parameter] = parallel_sampling(FBAmodel,parameter);

%estimate kinetic parameters in an ensemble
ensb = parallel_ensemble(FBAmodel,parameter,mc);

%solve ODE of model to steady state
sol = IntegrateModel(FBAmodel,ensb);
