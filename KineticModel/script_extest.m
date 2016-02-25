%script_extest
clc
addpath(genpath('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel'));
rxfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\glc_test.txt';
cnfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\glc_testC.txt';
%create model structure
[FBAmodel,parameter,variable,nrxn,nmetab] = modelgen(rxfname);

%obtain conentrations from file
[mc,FBAmodel,met] = readCNCfromFile(cnfname,FBAmodel);

% and/or

% sample initial metabolite concentrations for estimating kinetic parameters
rxn_add = {'GLCpts','NADH16'};
[mc,parameter,smp] = parallel_sampling(FBAmodel,parameter,'setupMetLP_glc',[],mc,rxn_add);

mc = iconcentration(FBAmodel,met,mc);
% Vup_struct.exO2 = 1000;%mmole/gDCW.h
% %fix flux uptakes 
% Vuptake = fixUptake(model,Vup_struct);
% if ~isfield(model,'Vuptake')
%     model.Vuptake = Vuptake;
% end

%get parameter estimates - estimate kinetic parameters in an ensemble
FBAmodel.rxn_add = rxn_add;

ensb = parallel_ensemble(FBAmodel,mc,parameter,rxn_add);

%serially solve ODE of model to steady state
if ensb{1,2}.feasible    
    sol = IntegrateModel(FBAmodel,ensb,ensb{1,1});
else
    error('No feasible model found');
end



