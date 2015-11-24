function ensb = parallel_ensemble(model,mc,pvec,rxn_add,rxn_excep,nmodels,smp)
if nargin<7
    smp={};
end

%# models from the ensemble
if nargin<6    
    nmodels=1;
end

%if no rxn exceptions for Vind in buildmodels
if nargin < 5
    rxn_excep = {};
end

if nargin < 4
    rxn_add = {};
end

%# metabolite concentration samples available
if ~isempty(smp)
    nsamples = size(smp,1);
else
    nsamples = 1;
end

if (nsamples==1 && nmodels==1)
    fprintf('\nGenerating a single model in the ensemble\n');
    ensb = cell(nmodels,2);
    model_ens = cell(nmodels,1);     
    ensb{1,1} = mc;
    model_ens{1} = buildmodels(model,pvec,mc,rxn_add,rxn_excep);    
elseif nsamples==1 && nmodels>1
    ensb = cell(nmodels,2);
    model_ens = cell(1,nmodels);
    %parallel capable
    parfor im = 1:nmodels
        ensb{im,1} = mc;
        model_ens{1,im} = buildmodels(model,pvec,mc,rxn_add,rxn_excep);
    end    
elseif nsamples>1 && nmodels>=1 %--------not supported yet
    ensb = cell(nsamples,nmodels,2);
    model_ens = cell(nsamples,nmodels);
    parfor ism = 1:nsamples        
        for im = 1:nmodels
            ensb{ism,im,1} = smp{ism,1};
            model_ens{ism,im} = buildmodels(model,smp{ism,2},smp{ism,1});
        end      
    end
end

if length(model_ens)>1
    for im = 1:length(model_ens)
        ensb{im,2} = model_ens{1,im};
    end
else
    ensb{1,2} = model_ens{1,1};
end

% if ndims(model_ens)>2
%     for ism = 1:nsamples
%         for im = 1:size(model_ens,1)
%             ensb{ism,im,2} = model_ens{ism,im};
%         end
%     end
% end
    



