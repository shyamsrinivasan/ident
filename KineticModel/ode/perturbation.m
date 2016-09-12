function [alloutsol,allpxeq,allpfeq] = perturbation(xeq,model,pvec,solverP,npts,type)
if nargin<6
    type = 'rnd';
end
if nargin<5
    npts = 1;
end

nmodels = size(pvec,2);nvar = size(xeq,1);

if nmodels>10
    parfor im = 1:nmodels
    end
else
    allpxeq = zeros(nvar,nmodels);
    allpfeq = zeros(length(model.rxns),nmodels);
    alloutsol(nmodels) = struct();
    for im = 1:nmodels
        % get npts perturbed points of xeq of type
        Nxeq = eqperturbation(xeq(:,im),npts,type);
        
        % solve each model strating from Nxeq
        [outsol,outss,xss,fss] =...
        solveAllpvec(model,pvec(im),Nxeq,solverP);
    
        % parse results for each perturbation for each models 
        alloutsol(im).t = outsol.t;
        alloutsol(im).y = outsol.y;
        alloutsol(im).flux = outsol.flux;
        allpxeq(:,im) = xss;
        allpfeq(:,im) = fss;
    end
end