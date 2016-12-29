% perturbation of points on the 2D manifold surface
function [] = Perturb2DstableManifold(man2Dpts,saddle,model,pvec,opts,tspanf,hfig,npts)
% 2Dmanpts - points on the 2D manifold represented as triangles from
% delaunay triangulation or from the known trajectories
nmanpts = size(man2Dpts,1);

% choose points in random that are away from each other
rndid = randi(nmanpts,npts,1);
startval = man2Dpts(rndid,:);

% perturb said points
posival = startval + repmat(1e-2*[0 0 1],npts,1);
negival = startval - repmat(1e-2*[0 0 1],npts,1);

% simulate from said points in forward time
for ival = 1:npts
    [xdynf1,xeq1] = solveODEonly(1,posival(ival,:)',model,pvec,opts,tspanf);
    [xdynf2,xeq2] = solveODEonly(1,negival(ival,:)',model,pvec,opts,tspanf);    
    
    % xdynr1 & xdynr2 plots
    FIGodetrajectories(xdynf1,posival(ival,:)',xeq1,2,[1 2 3],hfig);
    if xeq1(1)<xeq1(2)
        FIGmssEqIvalPerturbations(posival(ival,:)',xeq1,2,[1 2 3],hfig);
    end
    FIGodetrajectories(xdynf2,negival(ival,:)',xeq2,2,[1 2 3],hfig);
    if xeq2(1)<xeq2(2)
        FIGmssEqIvalPerturbations(negival(ival,:)',xeq2,2,[1 2 3],hfig);
    end
end

% check of final ss is the saddle node