function [model,Amat,Bmat,Cmat,Ns] = globalSensMatrix(model,nsampl,parScale)
[model.Coefficient,model.brate] = parameter_return(model.allpar,model);
%Coefficient
% [ngene,nreg] = size(model.RS);
[cr,cc] = find(model.Coefficient);
cf_ind = sub2ind(size(model.Coefficient),cr,cc);
ncoeff = length(cf_ind);
ceffind = 1:ncoeff;
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
nkcats = 1;
kcatind = 1+ninputs+ncoeff+nbr+ntr+nub+nb;

Ns = ncoeff+nbr+ntr+nub+nb+ninputs+nkcats;

pScale = zeros(Ns,2);
pScale(ceffind,:) = repmat(parScale(1,:),length(ceffind),1);
pScale(brind,[1,2]) = repmat(parScale(2,:),length(brind),1);
pScale(trind,[1,2]) = repmat(parScale(3,:),length(trind),1);
pScale(nbind,[1,2]) = repmat(parScale(4,:),length(nbind),1);
pScale(nubind,[1,2]) = repmat(parScale(5,:),length(nubind),1);
pScale(inpind(1),[1,2]) = parScale(6,:);
pScale(inpind(2),[1,2]) = parScale(7,:);
pScale(kcatind,[1,2]) = parScale(8,:);

pScale_lw = repmat(pScale(:,1),2,1)';
pScale_up = repmat(pScale(:,2),2,1)';
pScale_lw = repmat(pScale_lw,nsampl,1);
pScale_up = repmat(pScale_up,nsampl,1);

pd = makedist('Uniform');
% pSample = random(pd,nsampl,2*Ns);
parSampl = pScale_lw + (pScale_up - pScale_lw).*random(pd,nsampl,2*Ns);

Amat = parSampl(:,1:Ns);
Bmat = parSampl(:,Ns+1:2*Ns);
Cmat = zeros(nsampl,Ns,Ns);
for ic = 1:Ns    
    Cmat(:,:,ic) = Bmat;
    Cmat(:,ic,ic) = Amat(:,ic);
end
return