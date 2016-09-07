% kinetic ensemble modeling script any/all models
% testing with Kotte model of glucoeneogenesis
addpath(genpath('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel'));
% rxfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\data\gtoy1.txt';
% cnfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\data\gtoy1C.txt';
% rxfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\data\gtoy2.txt';
% cnfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\data\gtoy2C.txt';
rxfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\data\gtoy4.txt';
cnfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\data\gtoy4C.txt';
% create model structure
[model,parameter,variable,nrxn,nmetab] = modelgen(rxfname);

% obtain conentrations from file
[mc,model,met] = readCNCfromFile(cnfname,model);

% FBA or pFBA solution - optional
% fix flux uptakes for FBA solution
% Vup_struct.exAC = 20;%mmol/gDCW.h
Vupstruct.ACex = 20;
Vupstruct.ENZ1ex = 1;

% designate reactions for which uptake should not be zero in FBA
% ess_rxn = {'exH','exPI','exAC'};
essrxn = {'ACex'};

% assign initial fluxes and calculate FBA fluxes for direction
% FBAmodel.bmrxn = 14;
model.bmrxn = 3;
model = FBAfluxes(model,'fba',essrxn,Vupstruct,...
                     find(strcmpi(model.rxns,'bmex')));

rxnadd = {};
% Metabolite conecntrations 
% extracellular metabolites in M moles/L
% met.pi_e = 4e-2;

% sample initial metabolite concentrations for estimating kinetic parameters
[mc,parameter,smp] = parallel_sampling(model,parameter,'setupMetLP_gtoy',met,mc,rxnadd);
if isempty(mc)
    % multiple saples are being supplied
    mc = smp{1,1};
    parameter = smp{1,2};
end

% get parameter estimates - estimate kinetic parameters in an ensemble
rxnadd = {};
rxnexcep = {};
% FBAmodel.bmrxn = [];
ensb = parallel_ensemble(model,mc,parameter,rxnadd,rxnexcep,10);

% serially solve ODE of model to steady state
model.rxn_add = rxnadd;
model.rxn_excep = rxnexcep;

% setup model for integration 
[newmodel,newpvec,Nimc,solverP,flux,dXdt] =...
setupKineticODE(model,ensb,ensb{1,1},essrxn,Vupstruct,1000);

% solve only if models are feasible
if size(ensb,1)>1
    [outsol,outss,allxeq,allfeq] = solveAllpvec(newmodel,newpvec,Nimc,solverP);
else
    if ensb{1,2}.feasible    
        % integrate model
        [outsol,outss] = callODEsolver(newmodel,newpvec,Nimc,solverP);
    else
        error('No feasible model found');
    end
end

% time course plots
AllTimeCoursePlots(outsol,newmodel,{'pyr[c]','pep[c]','fdp[c]','ac[c]'},...
                                   {'ACt2r','FBP','PDHr','PYK'});
                               
% perturbations to steady states                              
