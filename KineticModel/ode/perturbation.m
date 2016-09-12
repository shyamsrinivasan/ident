function [alloutsol,allpxeq,allpfeq] =...
         perturbation(xeq,model,pvec,solverP,idx,npts,type)
if nargin<7
    type = 'rnd';
end
if nargin<6
    npts = 1;
end
if nargin<5
    idx = [];
else
    if iscell(idx)
        idx = cellfun(@(x)strcmpi(model.mets,x),idx,'UniformOutput',false);
        idx = cellfun(@(x)find(x),idx,'UniformOutput',false);
        idx = cell2mat(idx);
    end
end

nmodels = size(pvec,2);nvar = size(xeq,1);

if nmodels>10
    parfor im = 1:nmodels
    end
else
    allpxeq = zeros(nvar,nmodels);
    allpfeq = zeros(length(model.rxns),nmodels);
    alloutsol(nmodels) = struct();
    tstart = tic;
    fprintf('Perturbation of %d models...\n',nmodels);
    for im = 1:nmodels
        % get npts perturbed points of xeq of type
        Nxeq = eqperturbation(xeq(:,im),idx,npts,type);
        
        % solve each model strating from Nxeq
        [outsol,~,xss,fss] =...
        solveAllpvec(model,pvec(im),Nxeq,solverP);
    
        % parse results for each perturbation for each models 
        alloutsol(im).t = outsol.t;
        alloutsol(im).y = outsol.y;
        alloutsol(im).flux = outsol.flux;
        allpxeq(:,im) = xss;
        allpfeq(:,im) = fss;        
    end
    AllTimeCoursePlots(alloutsol,model,{'pyr[c]','pep[c]','fdp[c]','ac[c]'},...
                                   {'ACt2r','FBP','PDHr','PYK'});  
    fprintf('Time for complete perturbation simulations of %d model steady states: %4.3g\n',nmodels,toc(tstart));
end