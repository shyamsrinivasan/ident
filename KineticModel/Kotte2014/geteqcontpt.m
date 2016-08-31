% get points on equilibrium continuation curve close to prescribed values
% given as inputs
% s1 = varargin{1}; % output from MATCONT
% x1 = varargin{2}; % output from MATCONT
% pt = pvec(ap)
function [saddle,saddlepar,status] = geteqcontpt(s1,x1,pt,eps)
if nargin<4 || isempty(eps)
    eps = 1e-3;
end  
% get limit point boundaries
id = cat(1,s1.index);
if length(id)>2
    id = id(2:end-1);
end
% get continued variable values
eqcontvar = x1(end,id);
% pick smallest and largest parameter value and position
[mineqvar,minid] = min(eqcontvar);
[maxeqvar,maxid] = max(eqcontvar);
minid = id(minid);
maxid = id(maxid);
% check if input value is within bounds
status = 1;
if pt<mineqvar || pt>maxeqvar
    % if outside bounds - terminate and return a saddle node at the midpt
    saddle = [];
    saddlepar = [];
    while isempty(saddlepar)
        [saddle,saddlepar] = getsaddlenode(s1,x1,eps);
        eps = eps*10;
    end
    status = 0;
else
    % if pt within bounds - identify saddle closest to pt
    % get all data between minid and maxid
    reldata = x1(:,min([minid maxid]):max([minid maxid]));
    % search for parameter value closest to point 
    if any(abs(pt-reldata(end,:))<=eps)
        index = find(abs(pt-reldata(end,:))<=eps,1,'first');
        if ~isempty(index)
            saddle = reldata(:,index);
            saddlepar = saddle(end);
            saddle(end) = [];
        end
    else
        saddle = [];
        saddlepar = [];
        status = -1;
    end
end