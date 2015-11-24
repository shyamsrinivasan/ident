%script_toy
clc

addpath(genpath('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel'));
rxfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\N2mD_test.txt';

%create model structure
[FBAmodel,parameter,variable,nrxn,nmetab] = modelgen(rxfname);

%add rxn
% rxn.equation = 'atp[c] <==>';
% FBAmodel = addRxn(FBAmodel,parameter,rxn);

%fix flux uptakes for FBA solution
Vup_struct.exA = 20;%mmol/gDCW.h

%designate reactions for which uptake should be zero in FBA
ess_rxn = {'exA'};

%assign initial fluxes and calculate FBA fluxes for direction
% prxnid = find(strcmpi(FBAmodel.rxns,'exP'));
FBAmodel = FBAfluxes(FBAmodel,'',ess_rxn,Vup_struct);

%extracellular metabolites in M moles/L
% met.glc_e = 0.2;
met = struct([]);

%samwple initial metabolite concentrations for estimating kientic parameters
[mc,parameter,smp] = parallel_sampling(FBAmodel,parameter,met,'setupMetLP_toy',500);

%get one set of concentrations and coresponding delGr
% [x,delGr,assignFlag,vCorrectFlag] = getiConEstimate(FBAmodel);

%get more than one set of concentrations using ACHR sampling
% ACHRmetSampling(FBAmodel,1,500,200)

%Reactions to be considered cytosolic even if defined otherwise
%reactions to consider for kinetics other than Vind
rxn_add = {'Ain'};

%reactions not be considered cytosolic even if defined otherwise
rxn_excep = {};

FBAmodel.rxn_add = rxn_add;
FBAmodel.rxn_excep = rxn_excep;

%get parameter estimates - estimate kinetic parameters in an ensemble
ensb = parallel_ensemble(FBAmodel,mc,parameter,rxn_add,rxn_excep);

% x = initialsample(FBAmodel);

%change initial conditions to simulate a perturbation
% change_pos = [];

%solve ODE of model to steady state
if ensb{1,2}.feasible    
    sol = IntegrateModel(FBAmodel,ess_rxn,Vup_struct,ensb,ensb{1,1});
else
    error('No feasible model found');
end

