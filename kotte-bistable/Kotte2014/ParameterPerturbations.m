function varargout =...
ParameterPerturbations(model,pvec,xss,ivalpts,ivalid,xeqpts,eqid,ipt,...
                       tspan,colorSpec,opts,varargin)
if nargin < 9
    tspan = 1:2000;
end
if nargin < 10
    colorSpec = chooseColors(4,{'Green','Purple','Red','Orange'});
end
if nargin < 11
    opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
end
if nargin > 12
    switch nargin
        case 13
            hf1 = varargin{1};
            ha1 = [];
            f1 = 1;
            f2 = 0;
            f3 = 0;
        case 14
            hf1 = varargin{1};
            ha1 = varargin{2};
            f1 = 1;
            f2 = 0;
            f3 = 0;
        case 15
            hf1 = varargin{1};
            ha1 = varargin{2};
            hf2 = varargin{3};
            f1 = 1;
            f2 = 1;
            f3 = 0;
        case 16
            hf1 = varargin{1};
            ha1 = varargin{2};
            hf2 = varargin{3};
            ha2 = varargin{4};
            hf3 = [];
            f1 = 1;
            f2 = 1;
            f3 = 0;
        case 17
            hf1 = varargin{1};
            ha1 = varargin{2};
            hf2 = varargin{3};
            ha2 = varargin{4};
            hf3 = varargin{5};
            ha3 = [];
            f1 = 1;
            f2 = 1;
            f3 = 1;
        case 18
            hf1 = varargin{1};
            ha1 = varargin{2};
            hf2 = varargin{3};
            ha2 = varargin{4};
            hf3 = varargin{5};
            ha3 = varargin{6};
            f1 = 1;
            f2 = 1;
            f3 = 1;
    end    
else
    f1 = 0;
    f2 = 0;
    f3 = 0;
end


nvar = size(xss,1);
neq = size(xss,2);

ieq = 0;
while ieq < neq
    % perturbation from ieq ss
    ival1 = xss(:,ieq+1);
    [~,xeq1] =...
    solveODEonly(1,ival1,model,pvec,opts,tspan);    
    xeqpts(nvar*ieq+1:nvar*(ieq+1),ipt) = xeq1;
    ivalpts(nvar*ieq+1:nvar*(ieq+1),ipt) = ival1;
    if ival1(1)>ival1(2)
        ivalid(ieq+1,ipt) = 1;
    elseif ival1(2)>ival1(1)
        ivalid(ieq+1,ipt) = 2;    
    end
    
    if xeq1(1)>xeq1(2)
        % if pep > fdp - high ss                        
%         Point.MarkerFaceColor = colorSpec{1};   
%         Point.MarkerEdgeColor = colorSpec{1}; 
        eqid(ieq+1,ipt) = 1;
    elseif xeq1(2)>xeq1(1)
        % if pep < fdp - low ss                        
%         Point.MarkerFaceColor = colorSpec{2}; 
%         Point.MarkerEdgeColor = colorSpec{2};  
        eqid(ieq+1,ipt) = 2;
    end                    
%     Point.Marker = '.';
%     Point.MarkerSize = 25;
    % annotation
%     addanot.text = ['P' num2str(ipt)];
%     if f1
%         [hf1,ha1] =...
%         FIGmssEqIvalPerturbations(ival1,xeq1,2,[1 2],hf1,ha1,Point,addanot);
%     end
%     if f2
%         [hf2,ha2] =...
%         FIGmssEqIvalPerturbations(ival1,xeq1,2,[2 3],hf2,ha2,Point);
%     end
%     if f3
%         [hf3,ha3] =...
%         FIGmssEqIvalPerturbations(ival1,xeq1,2,[1 3],hf3,ha3,Point);
%     end
    ieq = ieq+1;
end

varargout{1} = ivalpts;
varargout{2} = ivalid;
varargout{3} = xeqpts;
varargout{4} = eqid;

if f1
    varargout{5} = hf1;
    varargout{6} = ha1;
end
if f2
    varargout{7} = hf2;
    varargout{8} = ha2;
end
if f3
    varargout{9} = hf3;
    varargout{10} = ha3;
end