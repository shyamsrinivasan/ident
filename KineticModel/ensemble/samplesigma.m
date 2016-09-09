function newK = samplesigma(model,M,K)
if nargin<3
    nmodels = 1;
end

S = model.S;
M = repmat(M,1,nmodels);
for irxn = 1:model.nt_rxn
    if irxn == Vind
        alls = S(:,irxn);alls(S(:,irxn)>0) = 0;
        allp = S(:,irxn);allp(S(:,irxn)<0) = 0;
        inKs = K(logical(alls),irxn);
        inKp = K(logical(allp),irxn);
        Kbkp = repmat([inKs;inKp],1,nmodels);
        ncmp = length(find(S(:,irxn)));
        sigma = random(makedist('Uniform',ncmp,nmodels));
        sats = sigma./(1-sigma);
        Kms = M(logical(alls),:)./...
        sats(1:length(find(alls)),:);
        Kmp = M(logical(allp),:)./...
        sats(length(find(alls)):length(find(alls))+length(find(allp)),:);
        Km = [Kms;Kmp];
        newK = Kbkp;
        if any(Kbkp(:,1))
            newK(Kbkp==1) = Km(Kbkp==1);
        end        
    end
end