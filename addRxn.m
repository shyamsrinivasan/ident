function model = addRxn(model,rxn)
if isfield(rxn,'equation')
    if iscell(rxn.equation)
        eqstr = rxn.equation;
    else
        eqstr{1} = rxn.equation;
    end
end
if isfield(rxn,'par')
    if iscell(rxn.par)
        allpar = rxn.par;
    else
        allpar{1} = rxn.par;
    end
else
    allpar = cell(length(eqstr),1);
end

for iadd = 1:length(eqstr)
%     ipt = 0;
    imetab = length(model.mets)+1;
    irxn = length(model.rxns)+1;

    [model,reverse] =...
    GetRxnInfo(model,imetab,irxn,allpar{iadd},eqstr{iadd});    
   
    model.reversible(irxn) = reverse;
    
    if isfield(rxn,'name')
        model.rxns{irxn} = rxn.name{iadd};
        model.enzs{irxn} = rxn.name{iadd};
    else
        model.rxns{irxn} = sprintf('newRxnID_%d',iadd);
        model.enzs{irxn} = sprintf('newRxnID_%d',iadd);
    end
    
    if isfield(rxn,'delSGr')
        model.delSGr(irxn) = rxn.delSGr(iadd);
    else
        model.delSGr(irxn) = 0;
    end
    
    if isfield(rxn,'delGlb')
    end
    
    if isfield(rxn,'delGub')
    end
    
    if isfield(rxn,'Keq')
    end
    
    
    
    [nmet,nrxn] = size(model.S);
    
    
    model.SI = [model.SI sparse(nmet,1)];
    
end
        