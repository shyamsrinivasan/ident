function jacobian = J_CKinetics(Vind)

nmet = size(S,1);
for irnx = 1:length(Vind)
    
    sbid = find(S(:,Vind(irxn))<0);
    prid = find(S(:,Vind(irxn))>0);
    
    Sb = -S(sbid,Vind(irxn));
    Sp = S(prid,Vind(irxn));
    
    Jnr = sparse(nmet,1);
    Jdr = sparse(nmet,1);
    
    for isb = 1:length(sbid)
        osid = setdiff(sbid,sbid(isb));
        %dNr/dS
        Jnr(sbid(isb)) = kfwd(Vind(irxn))*1/K(sbid(isb),Vind(irxn))*...
                         prod(mc(osid)./K(osid,Vind(irxn)));   
        %dDr/dS
        Kr = 1/K(sbid(isb),Vind(irxn));
        dr = Kr;
        for si = 2:Sb(sbid(isb))
            dr = dr+Kr*(si)*(1+mc(sbid(isb))./K(sbid(isb),Vind(irxn)))^(si-1);
        end
        dr_os = 1+mc(osid)./K(osid,Vind(irxn));
        for ios = 1:length(osid)
            for si = 2:Sb(osid(ios))
                dr_os(osid(ios)) = dr_os(osid(ios))+dr_os(osid(ios))^si;
            end
        end
        Jdr(sbid(isb)) = prod(dr_os)*dr;
    end
    
    for ipr = 1:length(prid)
        opid = setdiff(prid,prid(ipr));
        %dNr/dP
        Jnr(prid(ipr)) = -kbkw(Vind(irxn))*1/K(prid(ipr),Vind(irxn))*...
                         prod(mc(opid)./K(opid,Vind(irxn)));   
        %dDr/dP
        Kr = 1/K(prid(ipr),Vind(irxn));
        dr = Kr;
        for si = 2:Sp(prid(ipr))
            dr = dr+Kr*(si)*(1+mc(prid(ipr))./K(prid(ipr),Vind(irxn)))^(si-1);
        end
        dr_op = 1+mc(opid)./K(opid,Vind(irxn));
        for iop = 1:length(opid)
            for si = 2:Sp(opid(iop))
                dr_op(opid(iop)) = dr_op(opid(iop))+dr_op(opid(iop))^si;
            end
        end
        Jdr(prid(ipr)) = prod(dr_op)*dr;
    end
        
    if model.rev(Vind(irxn))
    elseif ~model.rev(Vind(irxn))
        
        Jnr = *
    end
    
end