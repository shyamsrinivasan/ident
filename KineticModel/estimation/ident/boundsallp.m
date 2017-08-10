function [newlb,newub] = boundsallp(lb,ub,data)

% km indices
kmid = [1,2,4,5,7,14,15];
% vmax/kcat indices
kcatid = [6,10,11,12,13];
% misc. concentrations
concid = 17;
% misc. parameters
miscid = [3,8,9,16];


newlb = zeros(data.nvar+1,1);
newub = zeros(data.nvar+1,1);

newlb(kmid) = .008;
newub(kmid) = 1;

newlb(kcatid) = .1;
newub(kcatid) = 2;

newlb(concid) = .1;
newub(concid) = .1;

newlb(miscid) = [1e6;1;.01;0];
newub(miscid) = [5e6;4;1;0];

newlb(data.idx) = [];
newub(data.idx) = [];

% if ~isempty(lb)
%     lb(data.varid) = newlb(data.varid);
%     newlb = lb;
% end
% if ~isempty(ub)
%     ub(data.varid) = newub(data.varid);
%     newub = ub;
% end






