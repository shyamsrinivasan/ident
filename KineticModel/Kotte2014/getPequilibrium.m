% get final equilibrium values for initial values perturbed from 
% unstable manifold
function getPequilibrium(ivals,model,pvec,opts,tspanf,hfig,id)

[npts,nvar] = size(ivals);

posp = ivals + repmat(1e-4*[1 1 1],npts,1);
negp = ivals - repmat(1e-4*[1 1 1],npts,1);

for jval = 1:npts
    [~,xeq1] = solveODEonly(1,posp(jval,:)',model,pvec,opts,tspanf);
    [~,xeq2] = solveODEonly(1,negp(jval,:)',model,pvec,opts,tspanf); 
    
    % equilibrium point plots
    if xeq1(1)<xeq1(2)
        FIGmssEqIvalPerturbations(posp(jval,:)',xeq1,2,id,hfig);
    end
    if xeq2(1)<xeq2(2)
        FIGmssEqIvalPerturbations(negp(jval,:)',xeq2,2,id,hfig);
    end
end
