% get points on equilibrium continuation curve close to prescribed values
% given as inputs
% s1 = varargin{1}; % output from MATCONT
% x1 = varargin{2}; % output from MATCONT
% pt = pvec(ap)
function [saddlepts,saddleparpts,statuspts] = geteqcontpt(s,pt,eps)
if nargin<3 || isempty(eps)
    eps = 1e-3;
end  

npts = size(fieldnames(s),1);
saddlepts = zeros(size(s.pt1.x1,1)-1,npts);
saddleparpts = zeros(1,npts);
statuspts = zeros(1,npts);
for ipt = 1:npts
    s1 = s.(['pt' num2str(ipt)]).s1;
    x1 = s.(['pt' num2str(ipt)]).x1;
    
    % get saddle point from here
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
    status = -1;
    if pt<mineqvar || pt>maxeqvar
        % if outside bounds - terminate and return a saddle node at the midpt        
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
        while status<0            
            if any(abs(pt-reldata(end,:))<=eps)
                index = find(abs(pt-reldata(end,:))<=eps,1,'first');
                if ~isempty(index)
                    saddle = reldata(:,index);
                    saddlepar = saddle(end);
                    saddle(end) = [];
                end
                status = 1;
            else
                saddle = [];
                saddlepar = [];                
                eps = eps*10;
            end
        end
    end
    % display(saddle); % debug
    if status>=0
        saddlepts(:,ipt) = saddle;
        saddleparpts(ipt) = saddlepar;
    end
    statuspts(ipt) = status;    
end % for




