function J = modelSparsity(model)
% sparsity pattern for J
[nx,~] = size(model.S);
J = sparse(0,nx);
% loop through metabolites
for ix = 1:nx
    % identify all reactions involved
    rxnid = find(model.S(ix,:)~=0);    
    % loop through reactions
    % identify all substrates, products and regulators
    [metid,~] = ind2sub([nx length(rxnid)],find(model.S(:,rxnid)~=0));
    [regid,~] = ind2sub([nx length(rxnid)],find(model.SI(:,rxnid)~=0));    
    Jin = sparse(1,unique([metid;regid]),1,1,nx);
    J = [J;Jin];
end
spy(J)