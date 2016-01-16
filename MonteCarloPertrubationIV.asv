function [sol] =...
MonteCarloPertrubationIV(model,ess_rxn,Vup_struct,ensb,initval,VMCpos,VMCneg,np)
%MC perturbation of initial values of the ODE for detecting multiplicity
if nargin<8
    np = 10;
    %number of perturbations
end
if nargin<7
    VMCneg = struct([]);
end
if nargin<6
    VMCpos = struct([]);
end
if nargin<5
    initval = zeros(model.nt_metab,1);
end
if nargin<4
    error('getiest:NoA','No parameter vector');
else
    pvec = ensb{1,2};
end
if nargin<3
    Vup_struct([]);
end
if nargin<2
    ess_rxn = {};
end
allMets = 1;%all metabolites

sol = cell(np,1);
finalSS = cell(np,1);

for ip = 1:np        
    if ~isempty(VMCpos)
        initval = changeInitialCondition(model,initval,VMCpos);
        allMets = 0;
    end
    if ~isempty(VMCneg)
        initval = changeInitialCondition(model,initval,[],VMCneg);
        allMets = 0;
    end
    if allMets
        initval = changeInitialCondition(model,initval,[],[],allMets);
    end    
    
    %initialize solver properties
    [model,solverP,saveData] = imodel(model,ess_rxn,Vup_struct,1e4);
    
    %calculate initial flux
    flux = iflux(model,pvec,initval.*model.imc);
    dXdt = ODEmodel(0,initval,[],model,pvec);
    
    %integrate model
    [sol{ip},finalSS{ip},status] = callODEsolver(model,pvec,initval,solverP);    
end

