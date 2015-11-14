function Vind = addToVind(model,rxn_add,rxn_excep)
Vind = model.Vind;

for irxn = 1:length(rxn_add)
    Vind = union(Vind,find(strcmpi(model.rxns,rxn_add{irxn})));
end

excep_ind = [];
for irxn = 1:length(rxn_excep)
    excep_ind = union(excep_ind,find(strcmpi(model.rxns,rxn_excep{irxn})));
end

Vind = setdiff(Vind,excep_ind);