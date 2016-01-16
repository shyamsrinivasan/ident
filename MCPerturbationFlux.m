function sol = MCPerturbationFlux(model,ess_rxn,Vup_struct,ensb,mc,initval,VMCpos,VMCneg,np)
%MC perturbation of fluxes of the ODE for detecting multiplicity
if nargin<9
    np = 10;
    %number of perturbations
end
if nargin<8
    VMCneg = struct([]);
end
if nargin<7
    VMCpos = struct([]);
end
if nargin<6
    initval = zeros(length(model.rxns),1);
end
if nargin<5
    mc = zeros(model.nt_metab,1);
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
allFlux = 1;%all fluxes

sol = cell(np,1);
finalSS = cell(np,1);

for ip = 1:np
    if ~isempty(VMCpos)
        initval = changeFlux(model,initval,VMCpos);
        allFlux = 0;
    end
    if ~isempty(VMCneg)
        initval = changeFlux(model,initval,[],VMCneg);
        allFlux = 0;
    end
    if allFlux
        initval = changeFlux(model,initval,[],[],allFlux);
    end    
    
    %initialize solver properties
    [model,solverP,saveData] = imodel(model,ess_rxn,Vup_struct,1e4);
    
    %calculate initial flux
    flux = iflux(model,pvec,mc.*model.imc);
    dXdt = ODEmodel(0,mc,[],model,pvec);
    
    %integrate model
    [sol{ip},finalSS{ip},status] = callODEsolver(model,pvec,initval,solverP);
end