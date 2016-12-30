function varargout = NumericalSeparatrix(varargin)
% numerical approximation of separatrix based on algorithm presented in
% Rudiger Seydel, Interdisciplinary Applied Mathematics, 1994
global model pvec colorSpec eps Axeq saddle
model = varargin{1};
pvec = varargin{2};
opts = varargin{3};
% 1. calculate equilibirum points
contvarid = varargin{4};
s1 = varargin{5};
x1 = varargin{6};
Axeq = varargin{7};
if nargin<8
    type = 'all';
else
    type = varargin{8};
end
if nargin<9
    tspanr = [0 -8.5];
else
    tspanr = varargin{9};
end
if nargin<10
    plotDim = 2;
else
    plotDim = varargin{10};
end
if nargin<11, 
    eps1 = [];
else
    eps1 = varargin{11};
end
ac = find(strcmpi(model.mets,'ac[e]'));

% 2. determine the saddle node/limit point
% calculate jacobian at each equilibirum point
eps = 1e-4;
% tspanr = [0,-18];
tspanf = 0:0.1:2000;
if plotDim == 2
    ht12f = [];
    ha12 = [];
    ht23f = [];
    ha23 = [];
    ht13f = [];
    ha13 = [];
    hf = [];
    ha = [];
    if nargout==4
        getdata = 1;
    end
end
if plotDim == 3
    hf3 = [];
    ha3 = [];
    hf = [];
    ha = [];
    if nargout==2
        getdata = 1;
    end
end
savexdynr = [];

% get saddle point
[saddle,saddlepar] = getsaddlenode(s1,x1,eps1);

% calculate characteristics at saddle point for confirmation
if ~isempty(saddle)
    unstablept = 0;
    pvec(contvarid) = saddlepar;
    model.PM(ac-length(saddle)) = saddlepar;
    [~,lambda,w] = getKotteJacobian(saddle,pvec,model);
    if any(real(lambda)>=0)
        unstablept = 1;
    end
else
    fprintf('No saddle points to work with\n');
    return
end

colorSpec = chooseColors(4,{'Green','Purple','Red','Navy','HotPink'});
Line.LineWidth = 2.0;
% select stable, unstable or all eigen vectors to be used
switch type
    case 'stable'
        tspanr1 = tspanr;
        xWs = calc1DWs(saddle,w,lambda,model,pvec,opts,tspanr1,tspanf,eps);
        Line.Color = colorSpec{1};   
        [hf,ha,hl1] = allTrajectoryfigures(plotDim,xWs,hf,ha,Line);
    case 'unstable'
        tspanr2 = tspanr;
        [xWus,xeq] = calc1DWus(saddle,w,lambda,model,pvec,opts,tspanr2,tspanf,eps);
        Line.Color = colorSpec{3};
        [hf,ha,hl2] = allTrajectoryfigures(plotDim,xWus,hf,ha,Line);       
        
        % separatrix perturbation close to saddle node
        Line1.Color = colorSpec{4};
        [hl4,hl5] = perturbSeparatrix(xWus,tspanf,opts,plotDim,hf,ha);                          
    otherwise
        fprintf('Using both stable and unstable eigen values for separatrix calculation\n');
%         wall = w;
%         lambdall = lambda;
%         w = [wall(:,real(lambda)>=0),wall(:,real(lambda)<0)];
%         lambda = [lambdall(real(lambdall)>=0);lambdall(real(lambdall)<0)];
        tspanr1 = tspanr(2,:);
        tspanr2 = tspanr(1,:);
        xWs = calc1DWs(saddle,w,lambda,model,pvec,opts,tspanr1,tspanf,eps);
        [xWus,xeq] = calc1DWus(saddle,w,lambda,model,pvec,opts,tspanr2,tspanf,eps);
        Line1 = Line;
        Line2 = Line;
        Line1.Color = colorSpec{1};   
        Line2.Color = colorSpec{3};
end         

% plot eq points on trajectories
if unstablept
    neqpts = size(Axeq,2);
    plotss = zeros(neqpts,1);
    Point.MarkerEdgeColor = [1 0 0];
    Point.MarkerFaceColor = [1 0 0];
    
    [hf,ha,hl3,plotss] = allPointfigures(plotDim,Axeq,xeq(:,1),hf,ha,Point,plotss);
    [hf,ha,hl] = allPointfigures(plotDim,Axeq,xeq(:,2),hf,ha,Point,plotss);
    
    nvar = length(saddle);     
    
    % stable ss vals from perturbing saddle node
%     [~,xeq] = solveODEonly(1,zval((iw-1)*nvar+1:iw*nvar,ival),model,pvec,opts,tspanf);
    % vector field around separatrix      

    if ~exist('hl1','var')
        hl1 = [];
    end
    if ~exist('hl2','var')
        hl2 = [];
    end

    if ~isempty(hl1)&&~isempty(hl2)
        legend([hl1 hl2 hl4 hl5 hl3(1) hl3(2)],'Separatrix(Stable Manifold)',...
               'Separatrix(Unstable Manifold)',...
               'Separatrix Perturbations',...
               'Separatrix Perturbations',...
               'Saddle Node','Steady States');
    elseif ~isempty(hl1)
        legend([hl1 hl4 hl5 hl3(1) hl3(2)],'Separatrix(Stable Manifold)',...               
               'Separatrix Perturbations',...
               'Separatrix Perturbations',...
               'Saddle Node','Steady States');
    elseif ~isempty(hl2)
        legend([hl2 hl4 hl5 hl3(1) hl3(2)],...
               'Separatrix(Unstable Manifold)',...
               'Separatrix Perturbations',...
               'Separatrix Perturbations',...
               'Saddle Node','Steady States');
    end        
%     legend('boxoff');
end

if plotDim == 2
    varargout{1} = hf(1);
    varargout{2} = hf(2);
    varargout{3} = hf(3);
    if getdata
        if strcmpi(type,'unstable')
            varargout{4} = xWus;
        elseif strcmpi(type,'stable')
            varargout{4} = xWs;
        end
    end
elseif plotDim == 3
    varargout{1} = hf;
    if getdata
        if strcmpi(type,'unstable')
            varargout{2} = xWus;
        elseif strcmpi(type,'stable')
            varargout{2} = xWs;
        end        
    end
end

function [hf,ha,hl,plotss] = allPointfigures(plotDim,Axeq,xeq,hf,ha,Point,plotss)
global saddle
ht12f = [];
ht23f = [];
ht13f = [];
ha12 = [];
ha23 = [];
ha13 = [];
neqpts = size(Axeq,2);

if any(any(abs(Axeq-repmat(real(xeq),1,neqpts))<=1e-8))
    [~,c] = find(abs(Axeq-repmat(real(xeq),1,neqpts))<=1e-8);
    c = unique(c);
    if ~plotss(c)                
        if plotDim == 2
            if ~isempty(hf)
                ht12f = hf(1);
                ht23f = hf(2);
                ht13f = hf(3);
            end
            if ~isempty(ha)            
                ha12 = ha(1);        
                ha23 = ha(2);        
                ha13 = ha(3);
            end
            % 2D stationary points (saddle/stable nodes)
            [~,~,hl] = FIGmssEqIvalPerturbations(saddle,real(xeq),...
                        2,[1 2],ht12f,ha12,Point);
            FIGmssEqIvalPerturbations(saddle,real(xeq),2,[1 3],...
                                      ht13f,ha13,Point);
            FIGmssEqIvalPerturbations(saddle,real(xeq),2,[2 3],...
                                      ht23f,ha23,Point);
        elseif plotDim == 3
            % 3D stationary points
            [~,~,hl] = FIGmssEqIvalPerturbations(saddle,real(xeq),...
                        2,[1 2 3],hf,ha,Point);
        end
        plotss(c) = 1;
    end
end

function [hf,ha,hl] = allTrajectoryfigures(plotDim,data,hf,ha,Line)
global saddle
ht12f = [];
ht23f = [];
ht13f = [];
ha12 = [];
ha23 = [];
ha13 = [];

if nargin<4
    Line = [];
end
% if ~isempty(ha)
if plotDim == 3
    [hf,ha,hl] =...
    FIGodetrajectories(real(data),saddle,saddle,2,[1 2 3],hf,ha,Line);
elseif plotDim == 2
    if ~isempty(hf)
        ht12f = hf(1);
        ht23f = hf(2);
        ht13f = hf(3);
    end
    if ~isempty(ha)            
        ha12 = ha(1);        
        ha23 = ha(2);        
        ha13 = ha(3);
    end
    [ht12f,ha12,hl1] =...
    FIGodetrajectories(real(data),saddle,saddle,2,[1 2],ht12f,ha12,Line);
    [ht23f,ha23] =...
    FIGodetrajectories(real(data),saddle,saddle,2,[2 3],ht23f,ha23,Line);
    [ht13f,ha13] =...
    FIGodetrajectories(real(data),saddle,saddle,2,[1 3],ht13f,ha13,Line);        
    hf = [ht12f;ht23f;ht13f];
    ha = [ha12;ha23;ha13];
end

function [hl1,hl2,xeq] = perturbSeparatrix(xdynr,tspanf,opts,plotDim,hf,ha)
global colorSpec eps model pvec Axeq saddle

if nargin<6, ha = []; end
if nargin<5, hf = []; end  
if nargin<4, plotDim = 2; end
hf12 = [];
hf23 = [];
hf13 = [];
ha12 = [];
ha23 = [];
ha13 = [];

Line.LineWidth = 2.0;
% perutbations from the separatrix - fwd integration only from mid point of
% separatrix curve
% if ~rem(size(xdynr,2)/2,2)
%     pid = size(xdynr,2)/2;
% else
%     pid = size(xdynr,2)/2+1/2;
% end
% get euclidean distance between saddle and xdynr points
edist = sqrt(sum((repmat(saddle,1,size(xdynr,2))-real(xdynr)).^2));
% select point within 5e-4 units from saddle
pid = find(edist<=5e-4,1,'last');

pival = xdynr(:,pid)+eps*[1;1;1]; % positive perturbation
[xdynf,xeq] = solveODEonly(1,pival,model,pvec,opts,tspanf);
if any(abs(Axeq(:,1)-real(xeq))<=1e-8)
    Line.Color = colorSpec{4};
elseif any(abs(Axeq(:,2)-real(xeq))<=1e-8)
    Line.Color = colorSpec{5};
end
if plotDim == 2
    if ~isempty(hf)
        hf12 = hf(1);
        hf23 = hf(2);
        hf13 = hf(3);
    end
    if ~isempty(ha)
        ha12 = ha(1);        
        ha23 = ha(2);        
        ha13 = ha(3);
    end
    [~,~,hl1] = FIGodetrajectories(real(xdynf),saddle,xeq,2,[1 2],hf12,ha12,Line);
    FIGodetrajectories(real(xdynf),saddle,xeq,2,[2 3],hf23,ha23,Line);
    FIGodetrajectories(real(xdynf),saddle,xeq,2,[1 3],hf13,ha13,Line);
elseif plotDim == 3
    if ~isempty(hf) && ~isempty(ha)
        [~,~,hl1] = FIGodetrajectories(real(xdynf),saddle,xeq,2,[1 2 3],hf,ha,Line);
    end
end

nival = xdynr(:,pid)-eps*[1;1;1]; % negative perturbation
[xdynf,xeq] = solveODEonly(1,nival,model,pvec,opts,tspanf);
if any(abs(Axeq(:,1)-real(xeq))<=1e-8)
    Line.Color = colorSpec{4};
elseif any(abs(Axeq(:,2)-real(xeq))<=1e-8)
    Line.Color = colorSpec{5};
end
if plotDim ==2
    if ~isempty(hf)
        hf12 = hf(1);
        hf23 = hf(2);
        hf13 = hf(3);
    end
    if ~isempty(ha)
        ha12 = ha(1);        
        ha23 = ha(2);        
        ha13 = ha(3);
    end
    [~,~,hl2] = FIGodetrajectories(real(xdynf),saddle,xeq,2,[1 2],hf12,ha12,Line);
    FIGodetrajectories(real(xdynf),saddle,xeq,2,[2 3],hf23,ha23,Line);
    FIGodetrajectories(real(xdynf),saddle,xeq,2,[1 3],hf13,ha13,Line);
elseif plotDim == 3
    if ~isempty(hf) && ~isempty(ha)
        [~,~,hl2] = FIGodetrajectories(real(xdynf),saddle,xeq,2,[1 2 3],hf,ha,Line);
    end
end


