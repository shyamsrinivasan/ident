%script_extest
% clc
addpath(genpath('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel'));
rxfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\ETC_test.txt';
cnfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\ETC_testC.txt';
% create model structure
[FBAmodel,parameter,variable,nrxn,nmetab] = modelgen(rxfname);

% obtain conentrations from file
[mc,FBAmodel,met] = readCNCfromFile(cnfname,FBAmodel);

% FBA or pFBA solution - optional
% fix flux uptakes for FBA solution
Vup_struct.exGLC = 20;%mmol/gDCW.h

% designate reactions for which uptake should not be zero in FBA
ess_rxn = {'exCO2','exH','exH2O','exPI','exO2','exGLC'};

% Optional - FBAmodel.Vss already has the requiste information from the
% excel file
% assign initial fluxes and calculate FBA fluxes for direction
FBAmodel.bmrxn = find(strcmpi(FBAmodel.rxns,'exPYR'));
FBAmodel = FBAfluxes(FBAmodel,'any',ess_rxn,Vup_struct);


% calculate delGr if concentrations cannot be sampled
rxn_add = {'GLCpts','NADH16','ATPS4r','CYTBD'};
% bounds = setupMetLP_g6p(FBAmodel,rxn_add,mc);
% [lnmc,assignFlag,delGr,vCorrectFlag] = assignConc(log(bounds.mc),FBAmodel,bounds);  
% parameter.delGr = delGr;

% and/or

% sample initial metabolite concentrations for estimating kinetic parameters
[mc,parameter,smp] = parallel_sampling(FBAmodel,parameter,'setupMetLP_red',met,mc,rxn_add);
if isempty(mc)
    % multiple saples are being supplied
    mc = smp{1,1};
    parameter = smp{1,2};
end

% [flux,vflux] = fluxATPS4r(FBAmodel,parameter,mc);

% mc = iconcentration(FBAmodel,met,exp(lnmc));
% Vup_struct.exO2 = 1000;%mmole/gDCW.h
% %fix flux uptakes 
% Vuptake = fixUptake(model,Vup_struct);
% if ~isfield(model,'Vuptake')
%     model.Vuptake = Vuptake;
% end

% get parameter estimates - estimate kinetic parameters in an ensemble
% FBAmodel.rxn_add = rxn_add;
rxn_add = {'GLCpts'};
rxn_excep = {'NADH16','ATPS4r','CYTBD','H2Ot'};
FBAmodel.bmrxn = [];
ensb = parallel_ensemble(FBAmodel,mc,parameter,rxn_add,rxn_excep);

%serially solve ODE of model to steady state
FBAmodel.rxn_add = rxn_add;
FBAmodel.rxn_excep = rxn_excep;
if ensb{1,2}.feasible    
    sol = IntegrateModel(FBAmodel,ensb,ensb{1,1});
else
    error('No feasible model found');
end



