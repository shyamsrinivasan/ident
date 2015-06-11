[nm,nr] = size(model.S);
[nm2,nr2] = size(FBAmodel.S);
newRxn = cell(nr2,1);
oldRxn = cell(nr2,1);
for ir = 1:nr
    for irxn = 1:nr2
        metid2 = FBAmodel.Metabolites(logical(FBAmodel.S(:,irxn)));
        metid = mets(logical(model.S(:,ir)));
        met_diff = setdiff(metid2,metid);
        if isempty(met_diff)
            newRxn{irxn} = model.rxns{ir};
            oldRxn{irxn} = FBAmodel.Enzyme{irxn};
        end
    end
end