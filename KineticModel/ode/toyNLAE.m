function dM = toyNLAE(mc,model,pvec)

nvar = size(mc,1);
dM = zeros(nvar,1);
PM = model.PM;
imc = model.imc;
% PM = cons(PM,mc);
% imc = cons(imc,mc);

% allmc = [mc.*imc;repmat(PM,1,size(mc,2))];
if ~isempty(PM)
    allmc = [mc.*imc;PM];
else
    allmc = mc.*imc;
end
% dM = cons(dM,allmc);

flux = iflux(model,pvec,allmc);

% metabolite consumption for biomass generation
fluxbm = BMKinetics(model,pvec,allmc,find(strcmpi(model.rxns,'biomass')));

% Cytosolic
dM(1:nvar) = (1./imc).*(model.S(1:nvar,:)*(flux)+fluxbm);
