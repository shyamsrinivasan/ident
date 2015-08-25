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
while  

%everything else