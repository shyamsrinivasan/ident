function dM = toyNLAE(mc,model,pvec)

nvar = size(mc,1);
dM = zeros(nvar,1);
PM = cons(model.Pimc,mc);
imc = cons(model.imc,mc);

allmc = [mc.*imc;repmat(PM,1,size(mc,2))];
dM = cons(dM,allmc);

flux = iflux(model,pvec,allmc);

% Cytosolic
dM(1:nvar) = (1./imc).*(model.S(1:nvar,:)*(flux*3600));
