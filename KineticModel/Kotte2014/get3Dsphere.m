function newivals = get3Dsphere(center,r,npts)
% get points around the centre within a sphere at a maximum distance r from it
% keep only positive points

theta = random('unif',0,2*pi,npts,1);
phi = random('unif',0,pi,npts,1);
r = random('unif',-r,r,npts,1);

x1 = r.*cos(theta).*sin(phi);
x2 = r.*sin(theta).*sin(phi);
x3 = r.*cos(phi);

newivals = repmat(center,npts,1) + [x1 x2 x3];

id = newivals(:,1)<0|newivals(:,2)<0|newivals(:,3)<0;
newivals(id,:) = [];
