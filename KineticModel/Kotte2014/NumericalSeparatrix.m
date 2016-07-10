function NumericalSeparatrix(varargin)
% numerical approximation of separatrix based on algorithm presented in
% Rudiger Seydel, Interdisciplinary Applied Mathematics, 1994
global model pvec colorSpec eps
model = varargin{1};
pvec = varargin{2};
opts = varargin{3};
% 1. calculate equilibirum points
eqpts = varargin{4};
eqcontvar = varargin{5};
contvarid = varargin{6};
s1 = varargin{7};
x1 = varargin{8};
ac = find(strcmpi(model.mets,'ac[e]'));

% 2. determine the saddle node/limit point
% calculate jacobian at each equilibirum point
saddle = [];
eps = 1e-4;
tspanr = [0,-15];
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
if any(abs(midpt-reldata(end,:))<=1e-3)
    index = find(abs(midpt-reldata(end,:))<=1e-3,1,'first');
    if ~isempty(index)
        saddle = reldata(:,index);
        saddlepar = saddle(end);
        saddle(end) = [];
    end
end

% calculate characteristics at saddle point for confirmation
if ~isempty(saddle)
    unstablept = 0;
    pvec(contvarid) = saddlepar;
    model.PM(ac-length(saddle)) = saddlepar;
    jacobian = getJacobian(saddle);
    [w,lambda] = eig(jacobian);
    lambda = diag(lambda);
    if any(real(lambda)>=0)
        unstablept = 1;
    end
else
    fprintf('No saddle points to work with\n');
    return
end

colorSpec = chooseColors(4,{'Green','Purple','Red','Orange'});
% reverse integration from saddle point in the direction of eigen vectors w
% at the saddle point
if unstablept
    p1execflag = 0;
    p2execflag = 0;
    Allxeq = [];
    for iw = 1:size(w,2)
        Point.MarkerEdgeColor = [1 0 0];
        Point.MarkerFaceColor = [1 0 0];
        % continuation parameter is already set
        zi = saddle+eps*w(:,iw);
        zj = saddle-eps*w(:,iw);
        xdynr = solveODEonly(1,zi,model,pvec,opts,tspanr);        
        [xdynf,xeq] = solveODEonly(1,zi,model,pvec,opts,tspanf);
        Line.Color = colorSpec{1};
        [hf3,ha3,hl1] = FIGodetrajectories(real(xdynr),saddle,saddle,2,[1 2 3],hf3,ha3,Line);
        [ht12fig,ha12] = FIGodetrajectories(real(xdynr),saddle,saddle,2,[1 2],ht12fig,ha12,Line);
        [ht23fig,ha23] = FIGodetrajectories(real(xdynr),saddle,saddle,2,[2 3],ht23fig,ha23,Line);
        [ht13fig,ha13] = FIGodetrajectories(real(xdynr),saddle,saddle,2,[1 3],ht13fig,ha13,Line);
        Line.Color = colorSpec{2};
        [hf3,ha3,hl2] = FIGodetrajectories(real(xdynf),saddle,saddle,2,[1 2 3],hf3,ha3,Line);
        [ht12fig,ha12] = FIGodetrajectories(real(xdynf),saddle,saddle,2,[1 2],ht12fig,ha12,Line);
        [ht23fig,ha23] = FIGodetrajectories(real(xdynf),saddle,saddle,2,[2 3],ht23fig,ha23,Line);
        [ht13fig,ha13] = FIGodetrajectories(real(xdynf),saddle,saddle,2,[1 3],ht13fig,ha13,Line);
        % plot equilibrium points
        if ~isempty(Allxeq)
            if ~any(abs(Allxeq-repmat(xeq,1,size(Allxeq,2)))<=1e-8)
                Allxeq = [Allxeq xeq];
                [~,~,hl3] = FIGmssEqIvalPerturbations(saddle,real(xeq),2,[1 2 3],hf3,ha3,Point);
                FIGmssEqIvalPerturbations(saddle,real(xeq),2,[1 2],ht12fig,ha12,Point);
                FIGmssEqIvalPerturbations(saddle,real(xeq),2,[1 3],ht13fig,ha13,Point);
                FIGmssEqIvalPerturbations(saddle,real(xeq),2,[2 3],ht23fig,ha23,Point);
            end
        else
            Allxeq = xeq;
            [~,~,hl3] = FIGmssEqIvalPerturbations(saddle,real(xeq),2,[1 2 3],hf3,ha3,Point);
            FIGmssEqIvalPerturbations(saddle,real(xeq),2,[1 2],ht12fig,ha12,Point);
            FIGmssEqIvalPerturbations(saddle,real(xeq),2,[1 3],ht13fig,ha13,Point);
            FIGmssEqIvalPerturbations(saddle,real(xeq),2,[2 3],ht23fig,ha23,Point);
        end
        % perturbations from separatrix
        if ~p1execflag
            hl4 = perturbSeparatrix(saddle,xdynr,tspanf,opts,hf3,ha3,...
                                    ht12fig,ha12,ht13fig,ha13,ht23fig,ha23);  
            p1execflag = 1;
        end
        
        xdynr = solveODEonly(1,zj,model,pvec,opts,[0 -9]);
        [xdynf,xeq] = solveODEonly(1,zj,model,pvec,opts,tspanf);
        Line.Color = colorSpec{1};
        [hf3,ha3] = FIGodetrajectories(real(xdynr),saddle,saddle,2,[1 2 3],hf3,ha3,Line);
        [ht12fig,ha12] = FIGodetrajectories(real(xdynr),saddle,saddle,2,[1 2],ht12fig,ha12,Line);
        [ht23fig,ha23] = FIGodetrajectories(real(xdynr),saddle,saddle,2,[2 3],ht23fig,ha23,Line);
        [ht13fig,ha13] = FIGodetrajectories(real(xdynr),saddle,saddle,2,[1 3],ht13fig,ha13,Line);
        Line.Color = colorSpec{2};
        [hf3,ha3] = FIGodetrajectories(real(xdynf),saddle,saddle,2,[1 2 3],hf3,ha3,Line);
        [ht12fig,ha12] = FIGodetrajectories(real(xdynf),saddle,saddle,2,[1 2],ht12fig,ha12,Line);
        [ht23fig,ha23] = FIGodetrajectories(real(xdynf),saddle,saddle,2,[2 3],ht23fig,ha23,Line);
        [ht13fig,ha13] = FIGodetrajectories(real(xdynf),saddle,saddle,2,[1 3],ht13fig,ha13,Line);
        % plot equilibrium points
        if ~isempty(Allxeq)
            if ~any(abs(Allxeq-repmat(xeq,1,size(Allxeq,2)))<=1e-8)
                Allxeq = [Allxeq xeq];
                FIGmssEqIvalPerturbations(saddle,real(xeq),2,[1 2 3],hf3,ha3,Point);
                FIGmssEqIvalPerturbations(saddle,real(xeq),2,[1 2],ht12fig,ha12,Point);
                FIGmssEqIvalPerturbations(saddle,real(xeq),2,[1 3],ht13fig,ha13,Point);
                FIGmssEqIvalPerturbations(saddle,real(xeq),2,[2 3],ht23fig,ha23,Point);
            end
        else
            Allxeq = xeq;
            FIGmssEqIvalPerturbations(saddle,real(xeq),2,[1 2 3],hf3,ha3,Point);
            FIGmssEqIvalPerturbations(saddle,real(xeq),2,[1 2],ht12fig,ha12,Point);
            FIGmssEqIvalPerturbations(saddle,real(xeq),2,[1 3],ht13fig,ha13,Point);
            FIGmssEqIvalPerturbations(saddle,real(xeq),2,[2 3],ht23fig,ha23,Point);
        end        
        % perturbations from separatrix
        if ~p2execflag
            perturbSeparatrix(saddle,xdynr,tspanf,opts,hf3,ha3,...
                              ht12fig,ha12,ht13fig,ha13,ht23fig,ha23);  
            p2execflag = 1;
        end
    end
    legend([hl1 hl2 hl4 hl3(1) hl3(2)],'Separatrix','Stable Manifold',...
           'Separatrix Perturbations','Saddle Node','Steady States');
end

function hl = perturbSeparatrix(saddle,xdynr,tspanf,opts,hf,ha,hf12,ha12,hf23,ha23,hf13,ha13)
% perutbations from the separatrix     
global colorSpec eps model pvec
if ~rem(size(xdynr,2)/2,2)
    pid = size(xdynr,2)/2;
else
    pid = size(xdynr,2)/2+1/2;
end
pival = xdynr(:,pid)+eps*[1;1;1];
[xdynf,xeq] = solveODEonly(1,pival,model,pvec,opts,tspanf);
Line.Color = colorSpec{4};
[~,~,hl] = FIGodetrajectories(real(xdynf),saddle,xeq,2,[1 2 3],hf,ha,Line);
FIGodetrajectories(real(xdynf),saddle,xeq,2,[1 2],hf12,ha12,Line);
FIGodetrajectories(real(xdynf),saddle,xeq,2,[2 3],hf23,ha23,Line);
FIGodetrajectories(real(xdynf),saddle,xeq,2,[1 3],hf13,ha13,Line);

nival = xdynr(:,pid)-eps*[1;1;1];
[xdynf,xeq] = solveODEonly(1,nival,model,pvec,opts,tspanf);
Line.Color = colorSpec{4};
FIGodetrajectories(real(xdynf),saddle,xeq,2,[1 2 3],hf,ha,Line);
FIGodetrajectories(real(xdynf),saddle,xeq,2,[1 2],hf12,ha12,Line);
FIGodetrajectories(real(xdynf),saddle,xeq,2,[2 3],hf23,ha23,Line);
FIGodetrajectories(real(xdynf),saddle,xeq,2,[1 3],hf13,ha13,Line);

function J = getJacobian(M)
global model pvec

% exact jacobian - for diagnostics only
Jxact = KottegivenJacobian(M,pvec,model);

% use ADMAT to calculate jacobians
admatfun = @(x)Kotte_givenNLAE(x,model,pvec);

xADMATobj = deriv(M,eye(3));
xADMATres = admatfun(xADMATobj);
% F = getval(xADMATres);
J = getydot(xADMATres); 

