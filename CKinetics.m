function flux = CKinetics(model,pvec,MC,Vind)
if nargin < 5
    useVmax = 1;
else
    useVmax = 0;
end
S = model.S;
nrxn = model.nt_rxn;
K = pvec.K;
kfwd = pvec.kcat_fwd;
kbkw = pvec.kcat_bkw;
Vmax = pvec.Vmax;

vflux = zeros(nrxn,1); 
flux = zeros(nrxn,1);
%eliminate consideration for excess cofators
%pi[c],pi[e],h[c],h[e],h2o[c]

pic = find(strcmpi(model.mets,'pi[c]'));
pie = find(strcmpi(model.mets,'pi[e]'));
hc = find(strcmpi(model.mets,'h[c]'));
he = find(strcmpi(model.mets,'h[e]'));
% q8 = find(strcmpi(model.mets,'q8[c]'));
% q8h2 = find(strcmpi(model.mets,'q8h2[c]'));
h2o = find(strcmpi(model.mets,'h2o[c]'));
co2 = find(strcmpi(model.mets,'co2[c]'));

for irxn = 1:length(Vind)
    sbid = S(:,Vind(irxn))<0;
    prid = S(:,Vind(irxn))>0;
    Sb = S(sbid,Vind(irxn));
    Sp = -S(prid,Vind(irxn));
    
    sbid([pic pie hc he h2o co2]) = 0;
    prid([pic pie hc he h2o co2]) = 0;
    
    if any(sbid) && any(prid)
        %Numerator - 1.6
        nr_flux = kfwd(irxn)*prod(MC(sbid)./K(sbid,Vind(irxn))) -...
                  kbkw(irxn)*prod(MC(prid)./K(prid,Vind(irxn)));
        %Denominator - 1.6
        %dr_sb
        dr_sb = 1+MC(sbid)./K(sbid,Vind(irxn));
        for j = 1:length(find(sbid))
            for si = 2:Sb(j)
                dr_sb(j) = dr_sb(j) + dr_sb(j)^si;                
            end
        end
        %dr_pr
        dr_pr = 1+MC(prid)./K(prid,Vind(irxn));
        for j = 1:length(find(prid))
            for si = 2:Sp(j)
                dr_pr(j) = dr_pr(j)+dr_pr(j)^si;
            end
        end
        dr_flux = prod(dr_sb)+prod(dr_pr)-1;
        vflux(Vind(irxn)) = nr_flux/dr_flux;
    else
        vflux(Vind(irxn)) = 0;
    end
    flux(Vind(irxn)) = Vmax(Vind(irxn))*vflux(Vind(irxn));
end
flux = flux(Vind);