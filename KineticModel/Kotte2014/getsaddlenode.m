function [saddle,saddlepar,index] = getsaddlenode(s1,x1,eps)
if nargin<3 || isempty(eps)
    eps = 1e-3;
end    
% s1 = varargin{1}; % output from MATCONT
% x1 = varargin{2}; % output from MATCONT
% get saddle point
id = cat(1,s1.index);
if length(id)>2
    id = id(2:end-1);
end
eqcontvar = x1(end,id);
% pick smallest and largest parameter value and position
[mineqvar,minid] = min(eqcontvar);
[maxeqvar,maxid] = max(eqcontvar);
minid = id(minid);
maxid = id(maxid);
% pick a midpoint parameter value
midpt = mineqvar+(maxeqvar-mineqvar)/2;
% get all data between minid and maxid
reldata = x1(:,min([minid maxid]):max([minid maxid]));
% search for parameter value closest to midpoint 
if any(abs(midpt-reldata(end,:))<=eps)
    index = find(abs(midpt-reldata(end,:))<=eps,1,'first');
    if ~isempty(index)
        saddle = reldata(:,index);
        saddlepar = saddle(end);
        saddle(end) = [];
    end
else
    saddle = [];
    saddlepar = [];
    index = [];
end