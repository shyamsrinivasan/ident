function flux = iflux(model,pvec,mc,flux)
if nargin<4
    flux = zeros(model.nt_rxn,1);
end
Vind = model.Vind;
Vex = model.Vex;
VFex = model.VFex;

%intracellular fluxes in Vind
flux(Vind) = CKinetics(model,pvec,mc,Vind);

%transport fluxes
flux(Vex) = TKinetics(model,pvec,mc,Vex);

%other fixed exchaged fluxes
flux(VFex) = EKinetics(model,mc,VFex);

%biomass
if mc(strcmpi(model.mets,'atp[c]'))>=1e-5
    flux(strcmpi(model.rxns,'atpm')) = 8.39;
else
    flux(strcmpi(model.rxns,'atpm')) = 0;
end

if flux(strcmpi(model.rxns,'atpm'))>=1e-5
    flux(model.bmrxn) = 0.1;
else
    flux(model.bmrxn) = 0;
end

