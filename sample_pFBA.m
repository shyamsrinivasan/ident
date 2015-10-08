function sample_pFBA(model)
%calcuate maximum growth rate
bounds.vl = zeros(model.nt_rxn,1);        
bounds.vl(logical(model.rev)) = -100;
bounds.vu = zeros(model.nt_rxn,1);          
bounds.vu(bounds.vu==0) = 100;
[vLPmax,vLPmin,model] = solveLP(model,bounds,model.bmrxn);
model.vl(model.c==1) = vLPmax.obj;

%call COBRA pFBA
model.lb = model.vl;
model.ub = model.vu;

load('C:\Users\shyam\Documents\E.Coli Metabolic Model_Feist_2007\EColi_iAF1260_Matlab\ecoli_core_model');
% model.c = model.c';
changeCobraSolver('ibm_cplex');
[GeneClasses RxnClasses modelIrrevFM] = pFBA(model, 'geneoption',0, 'tol',1e-7,'skipclass',1);
% modelIrrevFM.c = sparse(find(strcmpi('Biomass_Ecoli_core_w_GAM',modelIrrevFM.rxns)),1,1,size(modelIrrevFM.S,2),1);
x = optimizeCbModel(modelIrrevFM);
%reduce sum of all fluxes - pFBA

%convert to irreversible model
modelIrrev = convertIrreversible(model);

%add meatabolite to S
modelIrrev.S(end+1,:) = ones(size(modelIrrev.S(1,:)));
modelIrrev.b(end+1) = 0;
modelIrrev.mets{end+1} = 'fluxMeasure';

nm = size(modelIrrev.mets,1);

%add reaction to S matrix
modelIrrev.S(:,end+1) = sparse(find(strcmpi('fluxMeasure',modelIrrev.mets)),1,-1,nm,1);
modelIrrev.rxns{end+1} = 'netFlux';
modelIrrev.vl(end+1) = 0;
modelIrrev.vu(end+1) = Inf;
modelIrrev.rev(end+1) = 0;

prxnid = find(strcmpi('netFlux',modelIrrev.rxns));
bounds.vl = modelIrrev.vl;
bounds.vu = modelIrrev.vu;

[vLPmax,vLPmin,flag] = solveLP(modelIrrev,bounds,prxnid);
modelIrrev.vl(prxnid) = vLPmin.obj;
modelIrrev.vu(prxnid) = vLPmin.obj;
bounds.vl = modelIrrev.vl;
bounds.vu = modelIrrev.vu;

[vLPmax,vLPmin,flag] = solveLP(modelIrrev,bounds,model.bmrxn);


