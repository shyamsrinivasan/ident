function flux = EKinetics(model,mc,VFex)
flux = zeros(model.nt_rxn,1);

for irxn = 1:length(VFex)
    sbid = logical(model.S(:,VFex(irxn))<0);
    if ~model.Vuptake(VFex(irxn))
        flux(VFex(irxn)) = mc(sbid);
    else
        flux(VFex(irxn)) = -model.Vuptake(VFex(irxn));
    end
end
flux = flux(VFex);
