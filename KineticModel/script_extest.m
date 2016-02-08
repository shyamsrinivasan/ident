%script_extest
clc
addpath(genpath('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel'));
rxfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\o2_test.txt';

%create model structure
[FBAmodel,parameter,variable,nrxn,nmetab] = modelgen(rxfname);

ess_rxn = {'exCO2','exH','exH2O','exPI','exO2','exGLC'};

[model,bounds] = changebounds(model,ess_rxn);
ToyODEmodel(t,mc,data,model,pvec)