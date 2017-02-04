function perturbTrajectory(traj,model,pvec,opts,tspanf,hfig,id,npts)

% choose points in random that are away from each other
rndid = randi(nmanpts,npts,1);
startval = man2Dpts(rndid,:);