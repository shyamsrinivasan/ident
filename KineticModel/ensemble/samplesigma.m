function allK = samplesigma(model,M,K,kfwd,krev,Vind,nmodels)
if nargin<7
    nmodels = 1;
end

S = model.S;
Keq = model.Keq;
nrxn = model.nt_rxn;
vss = model.Vss;
rev = model.rev;

M = repmat(M,1,nmodels);
allK(nmodels) = struct();
[imet,jrxn] = find(K);
if nmodels>10
    parfor im = 1:nmodels
        allK(im).K = sparse(imet,jrxn,1,size(K,1),size(K,2));
        allK(im).kfwd = zeros(1,nrxn);
        allK(im).krev = zeros(1,nrxn);
    end
else
    for im = 1:nmodels
        allK(im).K = sparse(imet,jrxn,1,size(K,1),size(K,2));
        allK(im).kfwd = zeros(1,nrxn);
        allK(im).krev = zeros(1,nrxn);
    end
end

he = find(strcmpi(model.mets,'h[e]'));
hc = find(strcmpi(model.mets,'h[c]'));
h2o = find(strcmpi(model.mets,'h2o[c]'));

for irxn = 1:nrxn
    if ismember(irxn,Vind)
        % subsrates and products
        alls = S(:,irxn);alls(S(:,irxn)>0) = 0;
        allp = S(:,irxn);allp(S(:,irxn)<0) = 0;
        alls([hc he h2o],:) = 0;allp([hc he h2o],:) = 0;
        
        % substrate and product Kms backed up
        inKs = K(logical(alls),irxn);
        inKp = K(logical(allp),irxn);
        Kbkp = repmat([inKs;inKp],1,nmodels);
        ncmp = length(find(S(:,irxn)));
        
        % random estimation of nmodel sigma
        sigma = random(makedist('Uniform'),ncmp,nmodels);
        
        % saturation from sigma
        sats = sigma./(1-sigma);
        
        % Kms and Kmp from saturation and concentration
        Kms = M(logical(alls),:)./...
        sats(1:length(find(alls)),:);
        Kmp = M(logical(allp),:)./...
        sats(length(find(alls))+1:length(find(alls))+length(find(allp)),:);
        Km = [Kms;Kmp];
        
        % restore from backup
        newK = Kbkp;
        
        % only include Kms whose default value not specified in file
        if any(Kbkp(:,1))
            newK(Kbkp==1) = Km(Kbkp==1);
            % cant have zero Kms since M/0 = NaN
            newK(newK==0) = 1;
        end  
        
        % sample kcats for above newK                
        [newkf,newkr] = samplekcat(alls,allp,rev(irxn),kfwd(irxn),krev(irxn),newK,Keq(irxn));
        
        % assign values to parameter vectorsas fields in allK
        allK = assignKm(alls,allp,irxn,newK,newkf,newkr,allK,nmodels);   
    end
end

if nmodels>10
    parfor im = 1:nmodels 
        allK(im).kfwd = allK(im).kfwd';
        allK(im).krev = allK(im).krev';        
        allK(im).krev(~rev) = 0;
        allK(im).Vmax = ones(nrxn,1);
        allK(im).Vmax(vss==0) = 0;
    end
else
    for im = 1:nmodels    
        allK(im).kfwd = allK(im).kfwd';
        allK(im).krev = allK(im).krev';
        allK(im).krev(~rev) = 0;
        allK(im).Vmax = ones(nrxn,1);
        allK(im).Vmax(vss==0) = 0;
    end
end