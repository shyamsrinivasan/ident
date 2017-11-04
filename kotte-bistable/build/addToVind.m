function Vind = addToVind(rxnlst,Vind,rxn_add,rxn_excep)
if nargin<4
    rxn_excep = {};
end
% Vind = model.Vind;

for irxn = 1:length(rxn_add)
    Vind = union(Vind,find(strcmpi(rxnlst,rxn_add{irxn})));
end

excep_ind = [];
for irxn = 1:length(rxn_excep)
    excep_ind = union(excep_ind,find(strcmpi(rxnlst,rxn_excep{irxn})));
end

Vind = setdiff(Vind,excep_ind);
Vind = ToColumnVector(Vind);