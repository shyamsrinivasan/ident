function flux = EKinetics(model,flux,VFex)
% flux = zeros(model.nt_rxn,1);

for irxn = 1:length(VFex)
    sbid = logical(model.S(:,VFex(irxn))<0);
    vf = model.S(sbid,:);
    vdiff = setdiff(find(vf),VFex(irxn));
    if ~model.Vuptake(VFex(irxn))
        flux(VFex(irxn)) = model.S(sbid,vdiff)*flux(vdiff);
        %mc(sbid);
    else
        flux(VFex(irxn)) = -model.Vuptake(VFex(irxn));
    end
end
flux = flux(VFex);
