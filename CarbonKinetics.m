function [flux,vflux,vcup] = CarbonKinetics(model,pvec,mc,flux)
if nargin<4
    flux = zeros(model.nt_rxn,1);
end

vflux = zeros(model.nt_rxn,1);
Vmax = pvec.Vmax;

% h2o = find(strcmpi(model.mets,'h2o[c]'));
% pic = find(strcmpi(model.mets,'pi[c]'));
% pie = find(strcmpi(model.mets,'pi[e]'));
% hc = find(strcmpi(model.mets,'h[c]'));
% he = find(strcmpi(model.mets,'h[e]'));

vcup = find(strcmpi(model.rxns,'GLCpts'));

for iv = 1:length(vcup)
%     sbid = logical(model.S(:,Vex(irxn))<0);
%     prid = logical(model.S(:,Vex(irxn))>0);
%     
%     sbid([h2o pic pie hc he]) = 0;
%     prid([h2o pic pie hc he]) = 0;
    
    %vPTS
    [~,vflux(vcup(iv))] = CKinetics(model,pvec,mc,vcup(iv));
end

vflux(vcup) = scale_flux(vflux(vcup));
flux(vcup) = Vmax(vcup)*vflux(vcup);

% flux = flux(vcup);
vflux = vflux(vcup);