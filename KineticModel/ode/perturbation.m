function [outsol,outss,allxss,allfss] = perturbation(xeq,model,pvec,solverP,npts,type)
if nargin<6
    type = 'rnd';
end
if nargin<5
    npts = 1;
end

nmodels = size(pvec,2);

if nmodels>10
    parfor im = 1:nmodels
    end
else
    for im = 1:nmodels
        % get npts perturbed points of xeq of type
        Nxeq = eqperturbation(xeq(:,im),npts,type);
        
        % solve each model strating from Nxeq
        [outsol,outss,allxeq,allfeq,allJac,alllambda] =...
        solveAllpvec(model,pvec,Nxeq,solverP);
    
        % parse results for each perturbation for each models        
    end
end