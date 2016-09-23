function dM = toyNLAE(mc,model,pvec)

nvar = size(mc,1);
dM = zeros(nvar,1);
PM = model.PM;
imc = model.imc;
rho = model.rho;

Mint = model.Mint;
Mext = model.Mext;
Vind = model.Vind;
Vex = model.Vex;
biomass = strcmpi(model.mets,'biomass[e]');

% PM = cons(PM,mc);
% imc = cons(imc,mc);

% allmc = [mc.*imc;repmat(PM,1,size(mc,2))];
if ~isempty(PM)
    allmc = [mc.*imc;PM];
else
    allmc = mc.*imc;
end
% dM = cons(dM,allmc);

% flux in mmole/Lcw/h for Vind and Vex and mmole/Lc/h for VFex
[flux,fluxbm] = iflux(model,pvec,allmc);

% Cytosolic
dM(Mint) = (1./imc(Mint)).*(model.S(Mint,:)*(flux)+fluxbm(Mint));

% Extracellular
% change flux values for flux(Vind) and flux(Vex) mmole/Lcw/h -> mmole/Lc/h
flux([Vind Vex]) = changefluxunits(flux([Vind Vex]),allmc,model,1);

dM(Mext) = (1./imc(Mext)).*(model.S(Mext,:)*(flux));

% Biomass as in biomass[e] not in Mext
dM(biomass) = (1./imc(biomass)).*(model.S(biomass,:)*flux);
