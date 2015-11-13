%script_ecoli
clc

addpath(genpath('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel'));
rxfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\ecoliN1.txt';

%create model structure
[FBAmodel,parameter,variable,nrxn,nmetab] = modelgen(rxfname);

%add rxn
% rxn.equation = 'atp[c] <==>';
% FBAmodel = addRxn(FBAmodel,parameter,rxn);

%fix flux uptakes for FBA solution
Vup_struct.exGLC = 20;%mmol/gDCW.h
Vup_struct.exO2 = 1000;%mmole/gDCW.h

%designate reactions for which uptake should be zero in FBA
ess_rxn = {'exCO2','exH','exH2O','exPI','exO2','exGLC'};

%assign initial fluxes and calculate FBA fluxes for direction
FBAmodel = FBAfluxes(FBAmodel,'pfba',ess_rxn,Vup_struct);

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
if ensb{1,2}.feasible
    sol = IntegrateModel(FBAmodel,ensb,ensb{1,1});
else
    error('No feasible model found');
end
% sol = IntegrateModel(FBAmodel,ensb);
