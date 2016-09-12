function dM = toyNLAEemets(mc,model,pvec)
% NLAE for extracellular metabolites

nvar = size(mc,1);
nmets = model.nt_metab;
PM = model.PM;
imc = model.imc;
dM = zeros(nmets-nvar,1);

if ~isempty(PM)
    allmc = [mc.*imc;PM];
else
    allmc = mc.*imc;
end


flux = iflux(model,pvec,allmc);

% Extracellular
dM(1:(nmets-nvar)) = (1./imc).*(model.S(1:(nmets-nvar),:)*(flux));