function estimateVmax()
Vuptake = model.Vuptake;
%known fluxes
knindx = logical(Vuptake);%Vuptake only
n_known = length(find(knindx));
%initial mu guess
mu_g = 0.1;
%initial mu calculated
mu_c = biomass_flux(model,met,[],flux);
%unknown fluxes
while abs(mu_c-mu_g)<=1e-4
    %g6p - always first
    vglc = strcmpi(model.rxns,'exGLC');
    g6p = strcmpi(model.mets,'g6p[c]');
    vg6p = model.S(g6p,:)>0;
    [~,vcontr] = ConvenienceKinetics(model,pmeter,MC,vg6p);
    Vmax(vg6p) = flux(vglc)/vcontr;
    %rxnlist
    rxnlist = {''};
    
end

%everything else