function varargout = get2Dmanifoldpoints(points,model,pvec,tspanr,opts)

npoints = size(points,1);
% nvar = size(points,2);
% npts = length(tspanr)-2;
% allxdynr = zeros(nvar,(npts-2)*npoints);
allxdynr = [];
allx = zeros(npoints,length(tspanr));
ally = zeros(npoints,length(tspanr));
allz = zeros(npoints,length(tspanr));
for ip = 1:npoints
    xdynr = solveODEonly(1,points(ip,:)',model,pvec,opts,tspanr); 
    % find if any negative values
    r = find(xdynr(1,:)<0|xdynr(2,:)<0|xdynr(1,:)<0);
    if isempty(r)
        allx(ip,:) = xdynr(1,:);
        ally(ip,:) = xdynr(2,:);
        allz(ip,:) = xdynr(3,:);
    end
%     allxdynr(:,(ip-1)*npts+1:ip*npts) = xdynr(:,3:end);
%     allxdynr = [allxdynr xdynr(:,3:end)];
end

% x = allxdynr(1,:);
% y = allxdynr(2,:);
% z = allxdynr(3,:);

if nargout<=1
    varargout{1} = allxdynr;
elseif nargout>1
    varargout{1} = allx;
    varargout{2} = ally;
    varargout{3} = allz;
    varargout{4} = allxdynr;
end

