function model = createCOBRAmodel(oldmodel)
%split reversible reactions into 2 irreversible reactions
%Create COBRA model structure from model
rxns = size(oldmodel.S,2);
mets = size(oldmodel.S,1);
model.lb = zeros(rxns,1);
model.ub = zeros(rxns,1);
%Add lb and ub for fluxes
model.lb(model.rev == 0) = 0;
model.ub(model.rev == 0) = 1000;
model.lb(model.rev == 1) = -1000;
model.ub(model.rev == 1) = 1000;
%Add lb and ub for concentrations
model.xl
model.xu
%Add lb and ub for delG
model.dGl
model.dGu

model.c = sparse(model.bmrxn,1,1,rxns,1);
model.b = zeros(mets,1);
model.S = oldmodel.S;
return
    
    
% for irxn = rxns
%     rxnlength = size(model.S,2);
%     if strcmp(model.Rev{irxn},'rev')
%         model.S(:,rxnlength+1) = -model.S(:,irxn);
%         model.SI(:,rxnlength+1) = model.SI(:,irxn);
%         model.K(:,rxnlength+1) = model.K(:,irxn);
%         model.KIact(:,rxnlength+1) = model.KIact(:,irxn);
%         model.KIihb(:,rxnlength+1) = model.KIihb(:,irxn);
%         model.
%         model.rev(rxnlength+1) = 0;
%         %Other parameters
%     end    
% end