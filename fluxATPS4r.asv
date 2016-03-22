function vflux = fluxATPS4r(model,pvec,flux,mc)
% pmf based kinetics for ATPS4r

vatps4r = find(strcmpi(model.rxns,'ATPS4r'));
hc = find(strcmpi(model.mets,'h[c]'));
he = find(strcmpi(model.mets,'h[e]'));

sbid = model.S(:,Vind(irxn))<0;
prid = model.S(:,Vind(irxn))>0;   

sbid([hc he]) = 0;
prid([hc he]) = 0;

kfwd = pvec.kcat_fwd;
kbkw = pvec.kcat_bkw;
K = pvec.K;

vflux(vatps4r) = kfwd(vatps4r)

