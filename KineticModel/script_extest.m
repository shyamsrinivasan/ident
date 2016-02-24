%script_extest
clc
addpath(genpath('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel'));
rxfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\glc_test.txt';
cnfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\o2_testC.txt';
%create model structure
[FBAmodel,parameter,variable,nrxn,nmetab] = modelgen(rxfname);

%obtain conentrations from file
[mc,FBAmodel,met] = readCNCfromFile(cnfname,FBAmodel);

%setup problem for delGr calculation
bounds = setupMetLP_o2(FBAmodel,mc);

%delGr calculation
[lnmc,assignFlag,delGr,vCorrectFlag] = assignConc(log(bounds.mc),FBAmodel,bounds);
parameter.delGr = delGr;

%assign other metabolite concentrations from file

[mc,assignFlag] = iconcentration(FBAmodel,met,mc,assignFlag);

% Vup_struct.exO2 = 1000;%mmole/gDCW.h
% %fix flux uptakes 
% Vuptake = fixUptake(model,Vup_struct);
% if ~isfield(model,'Vuptake')
%     model.Vuptake = Vuptake;
% end

%get parameter estimates - estimate kinetic parameters in an ensemble
rxn_add = {'CYTBD','GLCpts','NADH16'};
FBAmodel.rxn_add = rxn_add;

ensb = parallel_ensemble(FBAmodel,mc,parameter,rxn_add);

%serially solve ODE of model to steady state
if ensb{1,2}.feasible    
    sol = IntegrateModel(FBAmodel,ensb,ensb{1,1});
else
    error('No feasible model found');
end



