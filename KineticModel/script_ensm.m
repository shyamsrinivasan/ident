% kinetic ensemble modeling script any/all models
% testing with Kotte model of glucoeneogenesis
addpath(genpath('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel'));
rxfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\gtoy1.txt';

%create model structure
[FBAmodel,parameter,variable,nrxn,nmetab] = modelgen(rxfname);

% FBA or pFBA solution - optional
% fix flux uptakes for FBA solution
Vup_struct.exAC = 20;%mmol/gDCW.h

% designate reactions for which uptake should not be zero in FBA
ess_rxn = {'exH','exPI','exAC'};

% assign initial fluxes and calculate FBA fluxes for direction
FBAmodel.bmrxn = 14;
FBAmodel = FBAfluxes(FBAmodel,'fba',ess_rxn,Vup_struct);

% Metabolite conecntrations 
% extracellular metabolites in M moles/L
met.pi_e = 4e-2;

% sample initial metabolite concentrations for estimating kinetic parameters
[mc,parameter,smp] = parallel_sampling(FBAmodel,parameter,'setupMetLP',met);
