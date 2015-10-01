function pvec = samplekcat(pvec,sbid,prid,irxn,mc)

delGr = pvec.delGr(irxn);     
kcatfwd = pvec.kcat_fwd(irxn);
kcatbkw = pvec.kcat_bkw(irxn);
K = pvec.K;

R = 0.008314; %kJ/mol.K
T = 298; %K
RT = R*T;

%normalized substrates
nrm_sb = prod(mc(sbid)./K(sbid,irxn));

%normalized product
nrm_pr = prod(mc(prid)./K(prid,irxn));

if isnan(kcatfwd)
    pvec.kcat_fwd(irxn) = nrm_pr/nrm_sb*kcatbkw*exp(-delGr/RT);
end
if isnan(kcatbkw)
    pvec.kcat_bkw(irxn) = nrm_sb/nrm_pr*kcatfwd*exp(delGr/RT);
end

pvec.kcat_fwd(irxn)*nrm_sb - pvec.kcat_bkw(irxn)*nrm_pr



    
    

