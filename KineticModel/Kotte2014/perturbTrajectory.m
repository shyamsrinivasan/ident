% choose points in random that are away from each other
function newivals = perturbTrajectory(traj,npts)
if nargin<2
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
skip = 0;
for jinterval = 1:ninterval
    if ~oldmin
        minx = min([traj(:,jinterval) traj(:,jinterval+1)],[],2);
    end
    maxx = max([traj(:,jinterval) traj(:,jinterval+1)],[],2);
    if any(abs(maxx-minx)>1e-4)
        newminx = min([minx maxx],[],2);
        newmaxx = max([minx maxx],[],2);        
    else
        oldmin = 1;
%         skip = skip+1;
%         ipts = ipts+ptspinterval;
        continue
    end
    % get numbers from a uniform distribution
    pd1 = makedist('Uniform','lower',newminx(1),'upper',newmaxx(1));
    pd2 = makedist('Uniform','lower',newminx(2),'upper',newmaxx(2));
    pd3 = makedist('Uniform','lower',newminx(3),'upper',newmaxx(3));
    r1 = random(pd1,ptspinterval,1);
    r2 = random(pd2,ptspinterval,1);
    r3 = random(pd3,ptspinterval,1);
%     if ipts == 
%         newivals(1:ptspinterval,:) = [r1 r2 r3];
    if ipts<newnpts
        newivals(ipts+1:ipts+ptspinterval,:) = [r1 r2 r3];
        ipts = ipts+ptspinterval;
    end    
end

newivals = newivals(1:ipts,:);


