function varargout = NumericalSeparatrix(varargin)
% numerical approximation of separatrix based on algorithm presented in
% Rudiger Seydel, Interdisciplinary Applied Mathematics, 1994
global model pvec colorSpec eps Axeq
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
    ht12fig = [];
    ha12 = [];
    ht23fig = [];
    ha23 = [];
    ht13fig = [];
    ha13 = [];
end
if plotDim == 3
    hf3 = [];
    ha3 = [];
end

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

% select stable, unstable or all eigen vectors to be used
switch type
    case 'stable'
        w = w(:,real(lambda)<0);
        lambda = lambda(real(lambda)<0);
        tspanr1 = tspanr;
        tspanr2 = [];
    case 'unstable'
        w = w(:,real(lambda)>=0);
        lambda = lambda(real(lambda)>=0);
        tspanr1 = [];
        tspanr2 = tspanr;
    otherwise
        fprintf('Using both stable and unstable eigen values for separatrix calculation\n');
        wall = w;
        lambdall = lambda;
        w = [wall(:,real(lambda)>=0),wall(:,real(lambda)<0)];
        lambda = [lambdall(real(lambdall)>=0);lambdall(real(lambdall)<0)];
        tspanr1 = tspanr(2,:);
        tspanr2 = tspanr(1,:);
end
    
colorSpec = chooseColors(4,{'Green','Purple','Red','Navy','HotPink'});
% Line.LineWidth = 2.5;
% reverse integration from saddle point in the direction of eigen vectors w
% at the saddle point
if unstablept
%     p1execflag = 0;
%     p2execflag = 0;
    neqpts = size(Axeq,2);
    plotss = zeros(neqpts,1);
    nvar = length(saddle);
    zval = zeros(nvar*size(w,2),2);
    for iw = 0:size(w,2)-1
        zval(iw*nvar+1:(iw+1)*nvar,:) =...
        [saddle+eps*w(:,iw+1) saddle-eps*w(:,iw+1)];
    end    
    Point.MarkerEdgeColor = [1 0 0];
    Point.MarkerFaceColor = [1 0 0];
    for ival = 1:size(zval,2);
        p1execflag = 0;
        for iw = 1:size(w,2)
            % separatrix curve
            if real(lambda(iw))<0
                if ~isempty(tspanr1)
                    xdynr1 =...
                    solveODEonly(1,zval((iw-1)*nvar+1:iw*nvar,ival),model,pvec,opts,tspanr1);
%                     xdynr1 = xdynr1(:,xdynr1(1,:)>=0);
%                     xdynr1 = xdynr1(:,xdynr1(2,:)>=0);
%                     xdynr1 = xdynr1(:,xdynr1(3,:)>=0);
                else
                    xdynr1 = [];
                end
            else
                xdynr1 = [];
            end
            if real(lambda(iw))>=0
                if ~isempty(tspanr2)
                    xdynr2 =...
                    solveODEonly(1,zval((iw-1)*nvar+1:iw*nvar,ival),model,pvec,opts,tspanr2);
%                     xdynr2 = xdynr2(:,xdynr2(1,:)>=0);
%                     xdynr2 = xdynr2(:,xdynr2(2,:)>=0);
%                     xdynr2 = xdynr2(:,xdynr2(3,:)>=0);
                else
                    xdynr2 = [];
                end
            else
                xdynr2 = [];
            end
            % stable ss vals
            [~,xeq] = solveODEonly(1,zval((iw-1)*nvar+1:iw*nvar,ival),model,pvec,opts,tspanf);
            % vector field around separatrix

            Line1.Color = colorSpec{1};   
            Line1.LineWidth = 2.0;
            Line2 = Line1;
            Line2.Color = colorSpec{3};
            % 2D trajectories (phase plane)
            if plotDim == 2
                % xdynr1
                if ~isempty(xdynr1)
                    [ht12fig,ha12,hl1] =...
                    FIGodetrajectories(real(xdynr1),saddle,saddle,2,[1 2],ht12fig,ha12,Line1);
                    [ht23fig,ha23] =...
                    FIGodetrajectories(real(xdynr1),saddle,saddle,2,[2 3],ht23fig,ha23,Line1);
                    [ht13fig,ha13] =...
                    FIGodetrajectories(real(xdynr1),saddle,saddle,2,[1 3],ht13fig,ha13,Line1);
%                 else
%                     hl1 = [];
                end
                % xdynr2
                if ~isempty(xdynr2)
                    [ht12fig,ha12,hl2] =...
                    FIGodetrajectories(real(xdynr2),saddle,saddle,2,[1 2],ht12fig,ha12,Line2);
                    [ht23fig,ha23] =...
                    FIGodetrajectories(real(xdynr2),saddle,saddle,2,[2 3],ht23fig,ha23,Line2);
                    [ht13fig,ha13] =...
                    FIGodetrajectories(real(xdynr2),saddle,saddle,2,[1 3],ht13fig,ha13,Line2);
%                 else
%                     hl2 = [];
                end
            end
            % 3D trajectories
            if plotDim == 3
                if ~isempty(xdynr1)
                    [hf3,ha3,hl1] =...
                    FIGodetrajectories(real(xdynr1),saddle,saddle,2,[1 2 3],hf3,ha3,Line1);
%                 else
%                     hl1 = [];
                end
                if ~isempty(xdynr2)
                    [hf3,ha3,hl2] =...
                    FIGodetrajectories(real(xdynr2),saddle,saddle,2,[1 2 3],hf3,ha3,Line2);
%                 else
%                     hl2 = [];
                end
            end
            % plot equilibrium points
            if any(any(abs(Axeq-repmat(real(xeq),1,neqpts))<=1e-8))
                [~,c] = find(abs(Axeq-repmat(real(xeq),1,neqpts))<=1e-8);
                c = unique(c);
                if ~plotss(c)                
                    if plotDim == 2
                        % 2D stationary points (saddle/stable nodes)
                        [~,~,hl3] = FIGmssEqIvalPerturbations(saddle,real(xeq),...
                                    2,[1 2],ht12fig,ha12,Point);
                        FIGmssEqIvalPerturbations(saddle,real(xeq),2,[1 3],...
                                                  ht13fig,ha13,Point);
                        FIGmssEqIvalPerturbations(saddle,real(xeq),2,[2 3],...
                                                  ht23fig,ha23,Point);
                    elseif plotDim == 3
                        % 3D stationary points
                        [~,~,hl3] = FIGmssEqIvalPerturbations(saddle,real(xeq),...
                                    2,[1 2 3],hf3,ha3,Point);
                    end
                    plotss(c) = 1;
                end
            end
            % perturbations from separatrix          
            if ~p1execflag
                Line1.Color = colorSpec{4};
                if plotDim == 2
                    if ~isempty(xdynr1)
                        hl4 = perturbSeparatrix(saddle,xdynr1,tspanf,opts,...
                                        ht12fig,ha12,ht23fig,ha23,ht13fig,ha13);  
                    end
                    if ~isempty(xdynr2)
                        hl4 = perturbSeparatrix(saddle,xdynr2,tspanf,opts,...
                                        ht12fig,ha12,ht23fig,ha23,ht13fig,ha13);  
                    end
                elseif plotDim == 3
                    if ~isempty(xdynr1)
                        hl4 = perturbSeparatrix(saddle,xdynr1,tspanf,opts,...
                                            [],[],[],[],[],[],plotDim,hf3,ha3);  
                    end
                    if ~isempty(xdynr2)
                        hl4 = perturbSeparatrix(saddle,xdynr2,tspanf,opts,...
                                        [],[],[],[],[],[],plotDim,hf3,ha3);  
                    end
                end
                p1execflag = 1;
            end
        end         
    end
    if ~exist('hl1','var')
        hl1 = [];
    end
    if ~exist('hl2','var')
        hl2 = [];
    end

    if ~isempty(hl1)&&~isempty(hl2)
        legend([hl1 hl2 hl4 hl3(1) hl3(2)],'Separatrix(Stable Manifold)',...
               'Separatrix(Unstable Manifold)',...
               'Separatrix Perturbations',...
               'Saddle Node','Steady States');
    elseif ~isempty(hl1)
        legend([hl1 hl4 hl3(1) hl3(2)],'Separatrix(Stable Manifold)',...               
               'Separatrix Perturbations',...
               'Saddle Node','Steady States');
    elseif ~isempty(hl2)
        legend([hl2 hl4 hl3(1) hl3(2)],...
               'Separatrix(Unstable Manifold)',...
               'Separatrix Perturbations',...
               'Saddle Node','Steady States');
    end        
%     legend('boxoff');
end

if plotDim == 2
    varargout{1} = ht12fig;
    varargout{2} = ht23fig;
    varargout{3} = ht13fig;
elseif plotDim == 3
    varargout{1} = hf3;
end

function [hl1,hl2,xeq] = perturbSeparatrix(saddle,xdynr,tspanf,opts,...
                                hf12,ha12,hf23,ha23,hf13,ha13,plotDim,hf,ha)
global colorSpec eps model pvec Axeq

if nargin<13, ha = []; end
if nargin<12, hf = []; end  
if nargin<11, plotDim = 2; end

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
    [~,~,hl2] = FIGodetrajectories(real(xdynf),saddle,xeq,2,[1 2],hf12,ha12,Line);
    FIGodetrajectories(real(xdynf),saddle,xeq,2,[2 3],hf23,ha23,Line);
    FIGodetrajectories(real(xdynf),saddle,xeq,2,[1 3],hf13,ha13,Line);
elseif plotDim == 3
    if ~isempty(hf) && ~isempty(ha)
        [~,~,hl2] = FIGodetrajectories(real(xdynf),saddle,xeq,2,[1 2 3],hf,ha,Line);
    end
end


