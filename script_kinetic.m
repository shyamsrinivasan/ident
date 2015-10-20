clc

addpath(genpath('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel'));
rxfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\ecoliN1.txt';

%create model structure
[FBAmodel,parameter,variable,nrxn,nmetab] = modelgen(rxfname);

%add rxn
% rxn.equation = 'atp[c] <==>';
% FBAmodel = addRxn(FBAmodel,parameter,rxn);

%assign initial fluxes and calculate FBA fluxes for direction
FBAmodel = FBAfluxes(FBAmodel,'pfba');

%samwple initial metabolite concentrations for estimating kientic parameters
[mc,parameter,smp] = parallel_sampling(FBAmodel,parameter,500);

%get one set of concentrations and coresponding delGr
% [x,delGr,assignFlag,vCorrectFlag] = getiConEstimate(FBAmodel);

%get more than one set of concentrations using ACHR sampling
% ACHRmetSampling(FBAmodel,1,500,200)

%get parameter estimates
ensb = parallel_ensemble(FBAmodel,mc,parameter);




% x = initialsample(FBAmodel);



%estimate kinetic parameters in an ensemble


%solve ODE of model to steady state
sol = IntegrateModel(FBAmodel,ensb);
