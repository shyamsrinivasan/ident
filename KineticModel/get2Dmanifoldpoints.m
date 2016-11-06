function allxdynr = get2Dmanifoldpoints(points,model,pvec,tspanr,opts)

npoints = size(points,1);
nvar = size(points,2);
npts = length(tspanr)-2;
allxdynr = zeros(nvar,(npts-2)*npoints);
for ip = 1:npoints
    xdynr = solveODEonly(1,points(ip,:)',model,pvec,opts,tspanr); 
    allxdynr(:,(ip-1)*npts+1:ip*npts) = xdynr(:,3:end);
end