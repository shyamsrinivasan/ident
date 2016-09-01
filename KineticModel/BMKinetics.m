function [flux,vflux] = BMKinetics(model,pvec,mc,Vid)
flux = zeros(model.nt_rxn,1);
vflux = zeros(model.nt_rxn,1);
Vmax = pvec.Vmax;
K = pvec.K;
KIact = pvec.KIact;
KIihb = pvec.KIihb;
S = model.S;
SI = model.SI;

he = find(strcmpi(model.mets,'h[e]'));
hc = find(strcmpi(model.mets,'h[c]'));
h2o = find(strcmpi(model.mets,'h2o[c]'));
fdp = strcmpi(model.mets,'fdp[c]');
pep = strcmpi(model.mets,'pep[c]');
% vectorize mc
vec_mc = repmat(mc,1,length(Vid));
% vec_mc = sparse(vec_mc);
vec_mc(~S(:,Vid)) = 0;

% vectorized mc for activation
actvec_mc = repmat(mc,1,length(Vid));
% actvec_mc = sparse(actvec_mc);
actvec_mc(SI(:,Vid)~=1) = 0;

% vectorized mc for inhibition
ihbvec_mc = repmat(mc,1,length(Vid));
% ihbvec_mc = sparse(ihbvec_mc);
ihbvec_mc(SI(:,Vid)~=-1) = 0;

% if rxn id is biomass reaction
if Vid == model.bmrxn
    ratio = 1+vec_mc(fdp)./K(fdp,Vid);
    vflux(Vid,:) = (ratio-1).*(ratio).^3./(ratio.^4+4e6.*...
                   (1+actvec_mc(pep)./KIact(pep,Vid)).^(-4));    
    vflux(Vid) = scale_flux(vflux(Vid,:));
    flux(Vid,:) = Vmax(Vid).*vflux(Vid,:);
end 

flux = flux(Vid,:);
vflux = vflux(Vid,:);