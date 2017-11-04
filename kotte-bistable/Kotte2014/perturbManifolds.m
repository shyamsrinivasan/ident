% perturbation of points on the 2D manifold surface
function [] = perturbManifolds(man2Dpts,saddle,model,pvec,opts,tspanf,hfig,id,npts)
% 2Dmanpts - points on the 2D manifold represented as triangles from
% delaunay triangulation or from the known trajectories
kid = man2Dpts(:,1)<0|man2Dpts(:,1)>5|man2Dpts(:,2)<0|...
      man2Dpts(:,2)>5|man2Dpts(:,3)<0|man2Dpts(:,3)>5;
  
man2Dpts(kid,:) = [];
nmanpts = size(man2Dpts,1);

% choose points in random that are away from each other
rndid = randi(nmanpts,npts,1);
startval = man2Dpts(rndid,:);

% perturb said points
posival = startval + repmat(1e-4*[1 1 1],npts,1);
negival = startval - repmat(1e-4*[1 1 1],npts,1);

% simulate from said points in forward time
for ival = 1:npts
    [xdynf1,xeq1] = solveODEonly(1,posival(ival,:)',model,pvec,opts,tspanf);
    [xdynf2,xeq2] = solveODEonly(1,negival(ival,:)',model,pvec,opts,tspanf);    
    
    % xdynr1 & xdynr2 plots
    FIGodetrajectories(xdynf1,posival(ival,:)',xeq1,2,id,hfig);
    if xeq1(1)<xeq1(2)
        FIGmssEqIvalPerturbations(posival(ival,:)',xeq1,2,id,hfig);
    end
    FIGodetrajectories(xdynf2,negival(ival,:)',xeq2,2,id,hfig);
    if xeq2(1)<xeq2(2)
        FIGmssEqIvalPerturbations(negival(ival,:)',xeq2,2,id,hfig);
    end
end

gcf
if length(id)==3
    axis([0 2.5 0 1.6 0 1.6]);
    view([116 22]);
elseif length(id) == 2   
    switch id(1)
        case 1
            xlim = [0 2.5];
        case 2
            xlim = [0 1.6];
        case 3
            xlim = [0 1.6];
    end
    switch id(2)
        case 1
            ylim = [0 2.5];
        case 2
            ylim = [0 1.6];
        case 3
            ylim = [0 1.6];
    end
    axis([xlim ylim]);    
end
% check of final ss is the saddle node