function NumericalSeparatrix(varargin)
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
if nargin<8, 
    eps1 = [];
else
    eps1 = varargin{8};
end
ac = find(strcmpi(model.mets,'ac[e]'));

% 2. determine the saddle node/limit point
% calculate jacobian at each equilibirum point
eps = 1e-4;
tspanr = [0,-8];
tspanf = 0:0.1:2000;
ht12fig = [];
ha12 = [];
ht23fig = [];
ha23 = [];
ht13fig = [];
ha13 = [];
hf3 = [];
ha3 = [];

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
% Line.LineWidth = 2.5;
% reverse integration from saddle point in the direction of eigen vectors w
% at the saddle point
if unstablept
    p1execflag = 0;
    p2execflag = 0;
    neqpts = size(Axeq,2);
    plotss = zeros(neqpts,1);
    Allxeq = [];
    for iw = 1:size(w,2)
        Point.MarkerEdgeColor = [1 0 0];
        Point.MarkerFaceColor = [1 0 0];
        % continuation parameter is already set
        zi = saddle+eps*w(:,iw);
        zj = saddle-eps*w(:,iw);
        % separatrix curve
        xdynr = solveODEonly(1,zi,model,pvec,opts,tspanr);   
        % stable manifold
        [~,xeq] = solveODEonly(1,zi,model,pvec,opts,tspanf);
        % vector field around separatrix
        
        Line.Color = colorSpec{1};   
        Line.LineWidth = 2.0;
        % 2D trajectories (phase plane)
        [ht12fig,ha12,hl1] = FIGodetrajectories(real(xdynr),saddle,saddle,2,[1 2],ht12fig,ha12,Line);
        [ht23fig,ha23] = FIGodetrajectories(real(xdynr),saddle,saddle,2,[2 3],ht23fig,ha23,Line);
        [ht13fig,ha13] = FIGodetrajectories(real(xdynr),saddle,saddle,2,[1 3],ht13fig,ha13,Line);        
%         Line.Color = colorSpec{2};        
%         [ht12fig,ha12,hl2] = FIGodetrajectories(real(xdynf),saddle,saddle,2,[1 2],ht12fig,ha12,Line);
%         [ht23fig,ha23] = FIGodetrajectories(real(xdynf),saddle,saddle,2,[2 3],ht23fig,ha23,Line);
%         [ht13fig,ha13] = FIGodetrajectories(real(xdynf),saddle,saddle,2,[1 3],ht13fig,ha13,Line);
        % plot equilibrium points
        if any(any(abs(Axeq-repmat(real(xeq),1,neqpts))<=1e-8))
            [~,c] = find(abs(Axeq-repmat(real(xeq),1,neqpts))<=1e-8);
            c = unique(c);
            if ~plotss(c)
                % 2D stationary points (saddle/stable nodes)
                [~,~,hl3] = FIGmssEqIvalPerturbations(saddle,real(xeq),...
                            2,[1 2],ht12fig,ha12,Point);
                FIGmssEqIvalPerturbations(saddle,real(xeq),2,[1 3],...
                                          ht13fig,ha13,Point);
                FIGmssEqIvalPerturbations(saddle,real(xeq),2,[2 3],...
                                          ht23fig,ha23,Point);
                plotss(c) = 1;
            end
        end
        % perturbations from separatrix
        if ~p1execflag
            Line.Color = colorSpec{4};
            hl4 = perturbSeparatrix(saddle,xdynr,tspanf,opts,...
                        ht12fig,ha12,ht23fig,ha23,ht13fig,ha13);  
            p1execflag = 1;
        end
        % separatrix curve
        xdynr = solveODEonly(1,zj,model,pvec,opts,[0 -8]);
        % stable manifold
        [~,xeq] = solveODEonly(1,zj,model,pvec,opts,tspanf);
        Line.Color = colorSpec{1};
        [ht12fig,ha12] = FIGodetrajectories(real(xdynr),saddle,saddle,2,[1 2],ht12fig,ha12,Line);
        [ht23fig,ha23] = FIGodetrajectories(real(xdynr),saddle,saddle,2,[2 3],ht23fig,ha23,Line);
        [ht13fig,ha13] = FIGodetrajectories(real(xdynr),saddle,saddle,2,[1 3],ht13fig,ha13,Line);
%         Line.Color = colorSpec{2};        
%         [ht12fig,ha12] = FIGodetrajectories(real(xdynf),saddle,saddle,2,[1 2],ht12fig,ha12,Line);
%         [ht23fig,ha23] = FIGodetrajectories(real(xdynf),saddle,saddle,2,[2 3],ht23fig,ha23,Line);
%         [ht13fig,ha13] = FIGodetrajectories(real(xdynf),saddle,saddle,2,[1 3],ht13fig,ha13,Line);
        % plot equilibrium points
        if any(any(abs(Axeq-repmat(real(xeq),1,neqpts))<=1e-8))
            [~,c] = find(abs(Axeq-repmat(real(xeq),1,neqpts))<=1e-8);
            c = unique(c);
            if ~plotss(c)
                % 2D stationary points (saddle/stable nodes)
                [~,~,hl3] = FIGmssEqIvalPerturbations(saddle,real(xeq),...
                            2,[1 2],ht12fig,ha12,Point);
                FIGmssEqIvalPerturbations(saddle,real(xeq),2,[1 3],...
                                          ht13fig,ha13,Point);
                FIGmssEqIvalPerturbations(saddle,real(xeq),2,[2 3],...
                                          ht23fig,ha23,Point);
                plotss(c) = 1;
            end
        end
        % perturbations from separatrix
        if ~p2execflag
            Line.Color = colorSpec{5};
            [~,hl5] = perturbSeparatrix(saddle,xdynr,tspanf,opts,...
                        ht12fig,ha12,ht23fig,ha23,ht13fig,ha13);  
            p2execflag = 1;
        end
    end
    legend([hl1 hl4 hl5 hl3(1) hl3(2)],'Separatrix',...
           'Separatrix Perturbations','Separatrix Perturbations',...
           'Saddle Node','Steady States');
    legend('boxoff');
end

function [hl1,hl2,xeq] = perturbSeparatrix(saddle,xdynr,tspanf,opts,...
                                hf12,ha12,hf23,ha23,hf13,ha13,hf,ha)
global colorSpec eps model pvec Axeq

if nargin<12, ha = []; end
if nargin<11, hf = []; end  

Line.LineWidth = 2.0;
% perutbations from the separatrix - fwd integration only
if ~rem(size(xdynr,2)/2,2)
    pid = size(xdynr,2)/2;
else
    pid = size(xdynr,2)/2+1/2;
end
pival = xdynr(:,pid)+eps*[1;1;1]; % positive perturbation
[xdynf,xeq] = solveODEonly(1,pival,model,pvec,opts,tspanf);
if any(abs(Axeq(:,1)-real(xeq))<=1e-8)
    Line.Color = colorSpec{4};
elseif any(abs(Axeq(:,2)-real(xeq))<=1e-8)
    Line.Color = colorSpec{5};
end
if ~isempty(hf) && ~isempty(ha)
    FIGodetrajectories(real(xdynf),saddle,xeq,2,[1 2 3],hf,ha,Line);
end
[~,~,hl1] = FIGodetrajectories(real(xdynf),saddle,xeq,2,[1 2],hf12,ha12,Line);
FIGodetrajectories(real(xdynf),saddle,xeq,2,[2 3],hf23,ha23,Line);
FIGodetrajectories(real(xdynf),saddle,xeq,2,[1 3],hf13,ha13,Line);

nival = xdynr(:,pid)-eps*[1;1;1]; % negative perturbation
[xdynf,xeq] = solveODEonly(1,nival,model,pvec,opts,tspanf);
if any(abs(Axeq(:,1)-real(xeq))<=1e-8)
    Line.Color = colorSpec{4};
elseif any(abs(Axeq(:,2)-real(xeq))<=1e-8)
    Line.Color = colorSpec{5};
end
if ~isempty(hf) && ~isempty(ha)
    FIGodetrajectories(real(xdynf),saddle,xeq,2,[1 2 3],hf,ha,Line);
end
[~,~,hl2] = FIGodetrajectories(real(xdynf),saddle,xeq,2,[1 2],hf12,ha12,Line);
FIGodetrajectories(real(xdynf),saddle,xeq,2,[2 3],hf23,ha23,Line);
FIGodetrajectories(real(xdynf),saddle,xeq,2,[1 3],hf13,ha13,Line);
