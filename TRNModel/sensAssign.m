function [model] = sensAssign(model,ASMP)
%Coefficient
[ngene,nreg] = size(model.RS);
[cr,cc] = find(model.Coefficient);
cf_ind = sub2ind(size(model.Coefficient),cr,cc);
ncoeff = length(cf_ind);
ceffind = 1:ncoeff;
% Coeff = sparse(ngene,nreg,1,cr,cc);
%Basal rate
nbr = size(model.RS,1);
brind = (1:nbr)+ncoeff;
%Translation
[tr,tc] = find(model.trate);
tr_ind = sub2ind(size(model.ptrate),tr,tc);
ntr = length(tr_ind);
trind = (1:ntr)+ncoeff+nbr;
% pTrate = sparse(nreg,ngene,1,tr,tc);
%binding & Unbinding
nub = length(find(model.Kub));
nubind = (1:nub)+ncoeff+nbr+ntr;
nb = length(find(model.Kb));
nbind = (1:nb)+ncoeff+nbr+ntr+nub;
%Inputs 
ninputs = 2;
inpind = (1:ninputs)+ncoeff+nbr+ntr+nub+nb;
%Kcat
% nkcats = 1;
kcatind = 1+ninputs+ncoeff+nbr+ntr+nub+nb;

model.Coefficient(cf_ind) = ASMP(ceffind);
model.brate(1:ngene) = ASMP(brind);
model.ptrate(tr_ind) = ASMP(trind);
model.Kb(model.Kb~=0) = ASMP(nbind);
model.Kub(model.Kub~=0) = ASMP(nubind);
model.gmax = ASMP(inpind(1));
model.Vuptake = ASMP(inpind(2));  
model.Kcat = zeros(model.nt_rxn,1);
model.Kcat(model.Kcat == 0) = ASMP(kcatind);
model.allpar = parameter_vector(model,length(model.Gene));
return           
    