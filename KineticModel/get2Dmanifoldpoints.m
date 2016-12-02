function varargout = get2Dmanifoldpoints(points,model,pvec,tspanr,opts)

npoints = size(points,1);
% nvar = size(points,2);
% npts = length(tspanr)-2;
% allxdynr = zeros(nvar,(npts-2)*npoints);
allxdynr = [];
for ip = 1:npoints
    xdynr = solveODEonly(1,points(ip,:)',model,pvec,opts,tspanr); 
%     allxdynr(:,(ip-1)*npts+1:ip*npts) = xdynr(:,3:end);
    allxdynr = [allxdynr xdynr(:,3:end)];
end

x = allxdynr(1,:);
y = allxdynr(2,:);
z = allxdynr(3,:);

chop_pos = 5;

r = find(x>chop_pos|x<0|y>chop_pos|y<0|z>chop_pos|z<0);
x(r) = [];
y(r) = [];
z(r) = [];

if nargout<=1
    varargout{1} = allxdynr;
elseif nargout>1
    varargout{1} = x;
    varargout{2} = y;
    varargout{3} = z;
end

