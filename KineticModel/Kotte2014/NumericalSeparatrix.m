function NumericalSeparatrix(varargin)
% numerical approximation of separatrix based on algorithm presented in
% Rudiger Seydel, Interdisciplinary Applied Mathematics, 1994
global model pvec
model = varargin{1};
pvec = varargin{2};
opts = varargin{3};
% 1. calculate equilibirum points
eqpts = varargin{4};
eqcontvar = varargin{5};
contvarid = varargin{6};
eigval = varargin{7};
ac = find(strcmpi(model.mets,'ac[e]'));

% 2. determine the saddle node/limit point
% calculate jacobian at each equilibirum point
saddle = [];
eps = 1e-4;
tspanr = [0,-10];
tspanf = 0:0.1:50;

% zi = eqpts(:,2)+eps*[1;0;0];
% xdynr = solveODEonly(1,zi,model,pvec,opts,tspan);
% zi = eqpts(:,2)+eps*[0;1;0];
% xdynr = solveODEonly(1,zi,model,pvec,opts,tspan);
% zi = eqpts(:,2)+eps*[0;0;1];
% xdynr = solveODEonly(1,zi,model,pvec,opts,tspan);
zi = eqpts(:,2)+eps*[1;1;1];
pvec(12) = eqcontvar(2);
model.PM(ac-length(zi)) = eqcontvar(2);
xdynr = solveODEonly(1,zi,model,pvec,opts,tspanr);

% zj = eqpts(:,2)-eps*[1;0;0];
% xdynr = solveODEonly(1,zj,model,pvec,opts,tspan);
% zj = eqpts(:,2)-eps*[0;1;0];
% xdynr = solveODEonly(1,zj,model,pvec,opts,tspan);
% zj = eqpts(:,2)-eps*[0;0;1];
% xdynr = solveODEonly(1,zj,model,pvec,opts,tspan);
zj = eqpts(:,2)-eps*[1;1;1];
xdynr = solveODEonly(1,zj,model,pvec,opts,tspanr);

% zi = eqpts(:,3)+eps*[1;0;0];
% xdynr = solveODEonly(1,zi,model,pvec,opts,tspan);
% zi = eqpts(:,3)+eps*[0;1;0];
% xdynr = solveODEonly(1,zi,model,pvec,opts,tspan);
% zi = eqpts(:,3)+eps*[0;0;1];
% xdynr = solveODEonly(1,zi,model,pvec,opts,tspan);
zi = eqpts(:,3)+eps*[1;1;1];
pvec(12) = eqcontvar(3);
model.PM(ac-length(zi)) = eqcontvar(3);
xdynr = solveODEonly(1,zi,model,pvec,opts,[0,-15]);

% zj = eqpts(:,3)-eps*[1;0;0];
% xdynr = solveODEonly(1,zj,model,pvec,opts,tspan);
% zj = eqpts(:,3)-eps*[0;1;0];
% xdynr = solveODEonly(1,zj,model,pvec,opts,tspan);
% zj = eqpts(:,3)-eps*[0;0;1];
% xdynr = solveODEonly(1,zj,model,pvec,opts,tspan);
zj = eqpts(:,3)-eps*[1;1;1];
xdynr = solveODEonly(1,zj,model,pvec,opts,[0,-15]);



for ieq = 1:size(eqpts,2)
    pvec(contvarid) = eqcontvar(ieq);
    jacobian = getJacobian(eqpts(:,ieq));
    % calcuate/obtain eigen value/vector 
    [w,lambda] = eig(jacobian);
    lambda = diag(lambda);
    % check for saddle node
    if any(real(lambda)>=0)
        saddle = eqpts(:,ieq);
    elseif ~isempty(saddle)
        saddle = [];    
    end
%     zi = saddle + eps*[1;0;0];
%     xdynr = solveODEonly(1,zi,model,pvec,opts,[0,-2]);
%     xdynf = solveODEonly(1,zi,model,pvec,opts,0:.1:50);
%     zj = saddle - eps*[1;0;0];
%     xdynr = solveODEonly(1,zj,model,pvec,opts,[0,-2]);
%     xdynf = solveODEonly(1,zj,model,pvec,opts,0:.1:50);
    
    % 3. get perturbed initial condition based on all eigen vectors
    if ~isempty(saddle)        
        for wi = 2:size(w,2)            
            zi = saddle + eps*w(:,wi);
            % 4a. integrate in reverse time to obtain the separatrix
            xdynr = solveODEonly(1,zi,model,pvec,opts,[0,-2]);
            % 4b. integration in forward time yields curves moving to the stable
            % equilibrium nodes
            xdynf = solveODEonly(1,zi,model,pvec,opts,0:.1:50);
            
            zj = saddle - eps*w(:,wi);
            % 4a. integrate in reverse time to obtain the separatrix
            xdynr = solveODEonly(1,zj,model,pvec,opts,[0,-2]);
            % 4b. integration in forward time yields curves moving to the stable
            % equilibrium nodes
            xdynf = solveODEonly(1,zj,model,pvec,opts,0:.1:50);
        end
    end
end

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

