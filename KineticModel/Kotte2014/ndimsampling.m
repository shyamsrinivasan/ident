function newivals = ndimsampling(centre,radius,npts)

pd = makedist('Uniform','lower',0,'upper',1);
r1 = random(pd,npts,3);
newivals = repmat(centre,npts,1) + r1.*repmat(radius,npts,1);