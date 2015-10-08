function flux = iflux(model,pvec,mc,flux)
if nargin<4
    flux = zeros(model.nt_rxn,1);
end

%intracellular fluxes in Vind
Vind = model.Vind;
flux(Vind) = CKinetics(model,pvec,mc,Vind);

%transport fluxes
Vex = model.Vex;
flux(Vex) = TKinetics(model,pvec,mc,Vex);

%other fixed exchaged fluxes
VFex = model.VFex;
flux(VFex) = EKinetics(model,flux,VFex);

%biomass
if mc(strcmpi(model.mets,'atp[c]'))>=1e-5
    flux(strcmpi(model.rxns,'atpm')) = 8.39;
else
    flux(strcmpi(model.rxns,'atpm')) = 0;
end

if flux(strcmpi(model.rxns,'atpm'))>=1e-5 &&...
    all(mc(logical(model.S(:,model.bmrxn)>0))>0)
    flux(model.bmrxn) = 0.1;
else
    flux(model.bmrxn) = 0;
end

