function delGr = delGsample(model)
if isfield('delGr',model)
    delGr = model.delGr;
else
    delGr = zeros(model.nt_rxn,1);
end

delGlb = model.delGlb;
delGub = model.delGub;
%sample delGr
for irxn = 1:model.nt_rxn
    %reaction direction
    %positive direction
    if delGlb(irxn)==0 && delGub(irxn)==0
        delGr(irxn) = 0;
    else
        if model.Vss(irxn) > 0
            delGr(irxn) = random(makedist('Uniform','lower',delGlb(irxn),'upper',0),1,1);
        elseif model.Vss(irxn) < 0
            delGr(irxn) = random(makedist('Uniform','lower',0,'upper',delGub(irxn)),1,1);
        elseif model.Vss(irxn) == 0
            delGr(irxn) = 0;        
        end
    end
end
