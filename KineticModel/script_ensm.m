% kinetic ensemble modeling script any/all models
% testing with Kotte model of glucoeneogenesis
addpath(genpath('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel'));
% rxfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\gtoy1.txt';
% cnfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\gtoy1C.txt';
% rxfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\gtoy2.txt';
% cnfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\gtoy2C.txt';
rxfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\gtoy3.txt';
cnfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\gtoy3C.txt';
% create model structure
[FBAmodel,parameter,variable,nrxn,nmetab] = modelgen(rxfname);

% obtain conentrations from file
[mc,FBAmodel,met] = readCNCfromFile(cnfname,FBAmodel);

% FBA or pFBA solution - optional
% fix flux uptakes for FBA solution
% Vup_struct.exAC = 20;%mmol/gDCW.h
Vupstruct.ACt2r = 20;
Vupstruct.ENZ1ex = 1;

% designate reactions for which uptake should not be zero in FBA
% ess_rxn = {'exH','exPI','exAC'};
essrxn = {'ACex'};

% assign initial fluxes and calculate FBA fluxes for direction
% FBAmodel.bmrxn = 14;
FBAmodel.bmrxn = 3;
FBAmodel = FBAfluxes(FBAmodel,'fba',essrxn,Vupstruct,...
                     find(strcmpi(FBAmodel.rxns,'bmex')));

rxnadd = {};
% Metabolite conecntrations 
% extracellular metabolites in M moles/L
% met.pi_e = 4e-2;

% sample initial metabolite concentrations for estimating kinetic parameters
[mc,parameter,smp] = parallel_sampling(FBAmodel,parameter,'setupMetLP_gtoy',met,mc,rxnadd);
if isempty(mc)
    % multiple saples are being supplied
    mc = smp{1,1};
    parameter = smp{1,2};
end

% get parameter estimates - estimate kinetic parameters in an ensemble
rxnadd = {};
rxnexcep = {};
% FBAmodel.bmrxn = [];
ensb = parallel_ensemble(FBAmodel,mc,parameter,rxnadd,rxnexcep);

% serially solve ODE of model to steady state
FBAmodel.rxn_add = rxnadd;
FBAmodel.rxn_excep = rxnexcep;
if ensb{1,2}.feasible    
    % setup model for integration 
    [newmodel,newpvec,Nimc,solverP,flux,dXdt] =...
    setupKineticODE(FBAmodel,ensb,ensb{1,1},essrxn,Vupstruct,15000);
    
    % integrate model
    [outsol,allxeq] = callODEsolver(newmodel,newpvec,Nimc,solverP);
%     sol = IntegrateModel(FBAmodel,ensb,ensb{1,1},ess_rxn,Vup_struct);
else
    error('No feasible model found');
end