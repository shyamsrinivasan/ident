function ensb = resample_ensemble(ensb,model,variable)
nmodels = length(fieldnames(ensb));
for isample = 1:nmodels
    mname = sprintf('model%d',isample);
    set = ensb.(mname);
    nrxn = size(set.K,2);
    for irxn = 1:nrxn
        rpind = model.S(:,irxn)~=0;
        rpK = set.K(rpind,irxn);
        rpC = variable.MC(rpind);
        %order of magnitude for Km based on C
        %concentration +/- 100
        if any(round(log10(rpK)-log10(rpC))>1)
%             respl = rpK(
            %resample these Ks
%             while log
        end
        
    end
end
return