function ensb = parallel_ensemble(model,pvec,mc)
 
if isstruct(mc)
    nsamples = length(fieldnames(mc));
else
    nsamples = 1;
end

smp = cell(nsamples,1);
if nsamples > 1
    %parallel possible
    for ism = 1:nsamples
        smp{ism} = buildmodels(model,pvec,mc);
    end
    
    for ism = 1:nsamples
        ensb.(sprintf('model%s',ism)) = smp{ism};
    end
else    
    ensb = buildmodels(model,pvec,mc);    
end