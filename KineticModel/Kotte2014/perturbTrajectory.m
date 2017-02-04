function perturbTrajectory(traj,model,pvec,opts,tspanf,hfig,id,npts)
if nargin<8
    npts = 5000;
end

% divide all npts points over the length of the trajectory 
% npts - number of points that need to be sampled ~10E3
[nvar,trajpts] = size(traj);
ninterval = trajpts-1;
ptspinterval = round(npts/ninterval);
newnpts = ptspinterval*ninterval;

newivals = zeros(newnpts,nvar);
ipts = 0;
oldmin = 0;
for jinterval = 1:ninterval
    if ~oldmin
        minx = min([traj(:,jinterval) traj(:,jinterval+1)],[],2);
    end
    maxx = max([traj(:,jinterval) traj(:,jinterval+1)],[],2);
    if any((maxx-minx)<1e-5)
        oldmin = 1;
        continue
    end
    % get numbers from a uniform distribution
    pd1 = makedist('Uniform','lower',minx(1),'upper',maxx(1));
    pd2 = makedist('Uniform','lower',minx(2),'upper',maxx(2));
    pd3 = makedist('Uniform','lower',minx(3),'upper',maxx(3));
    r1 = random(pd1,ptspinterval,1);
    r2 = random(pd2,ptspinterval,1);
    r3 = random(pd3,ptspinterval,1);
%     if ipts == 
%         newivals(1:ptspinterval,:) = [r1 r2 r3];
    if ipts>1 && ipts<newnpts
        newivals(ipts+1:ipts+ptspinterval,:) = [r1 r2 r3];
    end
    ipts = ipts+ptspinterval;
end

% choose points in random that are away from each other
