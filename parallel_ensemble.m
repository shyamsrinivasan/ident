function ensb = parallel_ensemble(model,mc,pvec,nmodels,smp)
if nargin<5
    smp={};
end
if nargin<4
    nmodels=1;
end
 
if ~isempty(smp)
    nsamples = size(smp,1);
else
    nsamples = 1;
end

if (nsamples==1 && nmodels==1)
    ensb = cell(nmodels,2);
    model_ens = cell(nmodels,1);     
    ensb{1,1} = mc;
    model_ens{1} = buildmodels(model,pvec,mc);    
elseif nsamples==1 && nmodels>1
    ensb = cell(nmodels,2);
    model_ens = cell(1,nmodels);
    %parallel capable
    parfor im = 1:nmodels
        ensb{im,1} = mc;
        model_ens{1,im} = buildmodels(model,pvec,mc);
    end    
elseif nsamples>1 && nmodels>=1
    ensb = cell(nsamples,nmodels,2);
    model_ens = cell(nsamples,nmodels);
    parfor ism = 1:nsamples        
        for im = 1:nmodels
            ensb{ism,im,1} = smp{ism,1};
            model_ens{ism,im} = buildmodels(model,smp{ism,2},smp{ism,1});
        end      
    end
end

if ndims(ensb)>2
    for ism = 1:nsamples
        for im = 1:size(model_ens,1)
            ensb{ism,im,2} = model_ens{ism,im};
        end
    end
else
    for im = 1:size(model_ens,1)
        ensb{im,2} = model_ens{1,im};
    end
end
    



