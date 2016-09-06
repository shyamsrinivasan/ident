function [outsol,outss,allxss,allfss] = solveAllpvec(model,pvec,ival,solverP)
% solve an ensemble of parameter vectors and post process results into a
% structure array

nmodels = size(pvec,2);
nvar = size(ival,1);
allsol = cell(nmodels,1);
allss = cell(nmodels,1);

fprintf('Initiating parallel solution to %d models...\n',nmodels);
tstart = tic;
parfor im = 1:nmodels
    if pvec(im).feasible
        [outsol,outss] = callODEsolver(model,pvec(im),ival,solverP);
        allsol{im} = outsol;
        allss{im} = outss;
    else
        allsol{im}.t = [];
        allsol{im}.y = [];
        allsol{im}.flux = [];
        allss{im}.t = [];
        allss{im}.y = [];
        allss{im}.flux = [];
    end
end
tstop = toc(tstart);
fprintf('\nParallel Integration time for %d models: %4.3g\n',nmodels,tstop);

outsol(nmodels) = struct();
outss(nmodels) = struct();
allxss = zeros(nvar,nmodels);
allfss = zeros(model.nt_rxn,nmodels);

for im = 1:nmodels
    if isfield(allsol{im},'t')
        outsol(im).t = allsol{im}.t;
    end
    if isfield(allsol{im},'y')
        outsol(im).y = allsol{im}.y;
    end
    if isfield(allsol{im},'flux')
        outsol(im).flux = allsol{im}.flux;
    end
    if isfield(allss{im},'t')
        outss(im).t = allss{im}.t;
    end
    if isfield(allss{im},'y')
        outss(im).y = allss{im}.y;
    end
    if isfield(allss{im},'flux')
        outss(im).flux = allss{im}.flux;
    end    
    % collect all ss together
    allxss(:,im) = outss(im).y;
    allfss(:,im) = outss(im).flux;
end
fprintf('Post processing overhead: %4.3g\n',toc(tstart)-tstop);