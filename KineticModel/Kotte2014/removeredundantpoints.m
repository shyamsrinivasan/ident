function [x,y,z] = removeredundantpoints(x,y,z,norm_tol)
% from Lyons 2014
if nargin<4
    norm_tol = 0.001;
end

% xbkp = x;
% ybkp = y;
% zbkp = z;

nxpts = size(x,2);
i = 1;
while i <= nxpts
    j = i+1;
    diffx = repmat(x(i),1,size(x(j:end),2))-x(j:end);
    diffy = repmat(y(i),1,size(y(j:end),2))-y(j:end);
    diffz = repmat(z(i),1,size(z(j:end),2))-z(j:end);
    
    eucdist = sqrt(diffx.^2+diffy.^2+diffz.^2);
    id = find(eucdist<norm_tol)+i;
    
    x(id) = [];
    y(id) = [];
    z(id) = [];
    i = i+1;
    
    nxpts = size(x,2);
end
%     
% x = xbkp;
% y = ybkp;
% z = zbkp;
% for i = 1:length(x)-1    
%     j = i + 1;
%     while  j <= length(x)
%         if sqrt((x(i)-x(j))^2 + (y(i)-y(j))^2 + (z(i)-z(j))^2) < norm_tol
%             x(j) = [];
%             y(j) = [];
%             z(j) = [];
%             j = j - 1;
%         end
%         j = j + 1;
%     end 
% end