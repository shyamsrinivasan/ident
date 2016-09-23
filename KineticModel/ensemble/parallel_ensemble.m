function [newp,mc] = parallel_ensemble(model,mc,pvec,rxnadd,rxnexcep,nmodels,smp)
if nargin<7
    smp={};
end
%# models from the ensemble
if nargin<6    
    nmodels=1;
end
%if no rxn exceptions for Vind in buildmodels
if nargin < 5
    rxnexcep = {};
end
if nargin < 4
    rxnadd = {};
end
%# metabolite concentration samples available
if ~isempty(smp)
    nsamples = size(smp,1);
else
    nsamples = 1;
end

% newK = samplesigma(model,mc,pvec.K,nmodels);

if nmodels>1
    tstart = tic;
    fprintf('\nGenerating %d parameter set in the ensemble...\n',nmodels);
    newp = buildmodels(model,pvec,mc,rxnadd,rxnexcep,nmodels);             
    fprintf('Parameter generation complete\n');
    fprintf('Time for generating %d model:%4.3g\n\n',nmodels,toc(tstart));   
else
    tstart = tic;
    fprintf('\nGenerating a single parameter set in the ensemble...\n');    
    newp = buildmodels(model,pvec,mc,rxnadd,rxnexcep,nmodels);    
    fprintf('Parameter generation complete\n');
    fprintf('Time for generating 1 model:%4.3g\n\n',toc(tstart));
end

% -------not supported
% if nsamples>1 && nmodels>=1 %--------not supported yet
%     ensb = cell(nsamples,nmodels,2);
%     model_ens = cell(nsamples,nmodels);
%     parfor ism = 1:nsamples        
%         for im = 1:nmodels
%             ensb{ism,im,1} = smp{ism,1};
%             model_ens{ism,im} = buildmodels(model,smp{ism,2},smp{ism,1});
%         end      
%     end
% end

% if length(model_ens)>1
%     for im = 1:length(model_ens)
%         ensb{im,2} = model_ens{1,im};
%     end
% else
%     ensb{1,2} = model_ens{1,1};
% end

% if ndims(model_ens)>2
%     for ism = 1:nsamples
%         for im = 1:size(model_ens,1)
%             ensb{ism,im,2} = model_ens{ism,im};
%         end
%     end
% end
    



