function [flux,vflux] = fluxATPS4r(model,pvec,mc,flux)
if nargin<4
    flux = zeros(length(model.rxns),1);
end
% pmf based kinetics for ATPS4r

vatps4r = find(strcmpi(model.rxns,'ATPS4r'));
hc = find(strcmpi(model.mets,'h[c]'));
he = find(strcmpi(model.mets,'h[e]'));

S = model.S;
Vmax = pvec.Vmax;

sbid = S(:,vatps4r)<0;
prid = S(:,vatps4r)>0;

sbid([hc he]) = 0;
prid([hc he]) = 0;

Sb = -S(sbid,vatps4r);
Sp = S(prid,vatps4r);

% kfwd = pvec.kcat_fwd;
% kbkw = pvec.kcat_bkw;
K = pvec.K(:,vatps4r);

% mche = 5*mc(hc)*exp(-(log(prod(mc(prid)))-log(prod(mc(sbid)))-log(model.Keq(vatps4r)))/4);


Katps = 220;
Kc = 2.5;
k1 = 5.13e3;
k2 = 2.16e3;

mc(hc) = 3.162e-9;
mc(he) = 1.585e-5;

delGr = log(prod(mc(prid)))-log(prod(mc(sbid)))-log(model.Keq(vatps4r))-...
        4*log(mc(hc)/mc(he));

x = (mc(he)./mc(hc)).*1/Katps;
D = 1 + x + x.^2 + x.^3 + x.^4;
kfwd = k2*1/(1+Kc)*1./D;
kbkw = k1*Kc/(1+Kc)*((x.^4)./D);
% kbkw = (kfwd/model.Keq(vatps4r))*prod(K(prid,1).^Sp)/prod(K(sbid,1).^Sb);
% kfwd = (kbkw*model.Keq(vatps4r))*prod(K(sbid,1).^Sb)/prod(K(prid,1).^Sp);
% 
% Vmax = 57;
% Km = 10^(-5.18);
% Ki = 7.34e-9;
% Kmp = Km*(1+mc(hc)*1e-3/Ki);
% nr_flux = Vmax*mc(he)*1e-3/(mc(he)*1e-3+Kmp)

nr_flux = kfwd*prod((mc(sbid)./K(sbid,1)).^Sb) -...
          kbkw*prod((mc(prid)./K(prid,1)).^Sp);
% ha = 1+d/ka;


if any(sbid) && any(prid)
    dr_sb = 1+mc(sbid)./K(sbid,1);
    for j = 1:length(find(sbid))
        for si = 2:Sb(j)
            dr_sb(j) = dr_sb(j) + dr_sb(j)^si;                
        end
    end
    dr_pr = 1+mc(prid)./K(prid,1);
    for j = 1:length(find(prid))
        for si = 2:Sp(j)
            dr_pr(j) = dr_pr(j)+dr_pr(j)^si;
        end
    end
else
    dr_sb = 0;
    dr_pr = 0;
end

if any(sbid) && any(prid)
    dr_flux = prod(dr_sb)+prod(dr_pr)-1;
    vflux(vatps4r) = scale_flux(nr_flux/dr_flux);
else
    vflux(vatps4r) = 0;
end
flux(vatps4r) = Vmax(vatps4r)*vflux(vatps4r);

flux = flux(vatps4r);
vflux = vflux(vatps4r);
end
     
   

