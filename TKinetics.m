function [flux,vflux] = TKinetics(model,pvec,mc,Vex)
flux = zeros(model.nt_rxn,1);
vflux = zeros(model.nt_rxn,1);
Vmax = pvec.Vmax;

h2o = find(strcmpi(model.mets,'h2o[c]'));
pic = find(strcmpi(model.mets,'pi[c]'));
pie = find(strcmpi(model.mets,'pi[e]'));
hc = find(strcmpi(model.mets,'h[c]'));
he = find(strcmpi(model.mets,'h[e]'));

vpts = find(strcmpi(model.rxns,'GLCpts'));
vthd2 = find(strcmpi(model.rxns,'THD2'));%nadph --->
vatps = find(strcmpi(model.rxns,'ATPS4r'));%atp ---> adp
vnad16 = find(strcmpi(model.rxns,'NADH16'));%nad ---> nadh
vcyt = find(strcmpi(model.rxns,'CYTBD'));%o2 --->
vspl = [vpts vthd2 vatps vnad16 vcyt];

nadph = strcmpi(model.mets,'nadph[c]');
nadp = strcmpi(model.mets,'nadp[c]');
nad = strcmpi(model.mets,'nad[c]');
nadh = strcmpi(model.mets,'nadh[c]');
atp = strcmpi(model.mets,'atp[c]');
adp = strcmpi(model.mets,'adp[c]');
o2 = strcmpi(model.mets,'o2[c]');

for irxn = 1:length(Vex)
    sbid = logical(model.S(:,Vex(irxn))<0);
    prid = logical(model.S(:,Vex(irxn))>0);
    
    %eliminate water from all reactions
    if ~strcmpi(model.rxns{Vex(irxn)},'h2ot')
        sbid(h2o) = 0;
        prid(h2o) = 0;
    end
        
    sbid([pic pie hc he]) = 0;
    prid([pic pie hc he]) = 0;
    
    if ~isnan(pvec.kcat_fwd(Vex(irxn)))
        kfwd = pvec.kcat_fwd(Vex(irxn));
    else
        kfwd = 1;
    end
    if ~isnan(pvec.kcat_bkw(Vex(irxn)))
        kbkw = pvec.kcat_bkw(Vex(irxn));
    else
        kbkw = 1;
    end
    
    if ~ismember(Vex(irxn),vspl)
        if model.rev(Vex(irxn))        
            vflux(Vex(irxn)) = kfwd*prod(mc(sbid)) - kbkw*prod(mc(prid));
        elseif ~model.rev(Vex(irxn))
            vflux(Vex(irxn)) = kfwd*prod(mc(sbid));
        end   
    end
    
    %other fluxes
    if Vex(irxn)==vpts
        [~,vflux(Vex(irxn))] = CKinetics(model,pvec,mc,Vex(irxn));
    end
    if Vex(irxn)==vthd2
        vflux(Vex(irxn)) = CKinetics(model,pvec,mc,Vex(irxn));
        %kfwd*mc(nadph)-kbkw*mc(nadp);
    end
    if Vex(irxn)==vatps
        vflux(Vex(irxn)) = CKinetics(model,pvec,mc,Vex(irxn));
        %kfwd*mc(atp)-kbkw*mc(adp);
    end
    if Vex(irxn)==vnad16
        vflux(Vex(irxn)) = CKinetics(model,pvec,mc,Vex(irxn));
        %kfwd*mc(nad)-kbkw*mc(nadh);
    end
    if Vex(irxn)==vcyt
        vflux(Vex(irxn)) = kfwd*(mc(o2))^0.5;
    end   
    vflux(Vex(irxn)) = scale_flux(vflux(Vex(irxn)));
    flux(Vex(irxn)) = Vmax(Vex(irxn))*vflux(Vex(irxn));
end
flux = flux(Vex);
vflux = vflux(Vex);
