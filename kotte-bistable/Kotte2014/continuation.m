function [contpt,varargout] = continuation(model,pvec,lambda,ap,iguess,opts)
if nargin<6
    opts = optimoptions(@fsolve,'TolX',1e-6,'TolFun',1e-6,'Display','iter');
end

% npts = size(iguess,1);
% contpt = zeros(npts,size(iguess,2));
% 
% for i = 1:npts
%     if sum(isnan(iguess(i,:)))==0
%         pvec(ap) = 
%         sysfun = @(x)Kotte_givenNLAE(x,model,pvec);
%         [fpsol,fval,flag] = fsolve(sysfun,iguess(i,:),opts);
%         fval = max(fval);
%     else
%         fpsol = NaN;
%         fval = 1;
%     end
%     contpt(i,:) = fpsol;
%     if exitflag~=1 || abs(fval)>opts
%         contpt(i,:) = NaN;
%     end
% end

ac = find(strcmpi(model.mets,'ac[e]'));
nvar = size(iguess,1);
% npts = size(iguess,1);
npar = length(lambda);
alleqsol = zeros(nvar,2*npar);
alleqflux = zeros(5,2*npar);
% eqpttype = zeros(npts,npar);
alleigval = zeros(nvar,2*npar);
% continuation algorithm
% set parameter value loop
for ip = 1:npar 
    model.PM(ac-nvar) = lambda(ip);
    pvec(ap) = lambda(ip);
    [~,neweqpt,~,eqflux] = solveODEonly(1,iguess,model,pvec,opts,1:2000);

    % store neweqpt
    alleqsol(:,ip) = neweqpt;
    alleqflux(:,ip) = eqflux;
    
    % get stability info
    [~,eigval] = KotteStabilityInfo(neweqpt',model,pvec);
    alleigval(:,ip) = eigval;
    
    % set initial value loop (from different solution branches ?)
    iguess = alleqsol(:,ip);
end

ipt = npar;
for ip = npar:-1:1
    ipt = ipt+1;
    model.PM(ac-nvar) = lambda(ip);
    pvec(ap) = lambda(ip);
    [~,neweqpt,~,eqflux] = solveODEonly(1,iguess,model,pvec,opts,1:2000);
    
    % store neweqpt
    alleqsol(:,ipt) = neweqpt;
    alleqflux(:,ipt) = eqflux;
    
    % get stability info
    [~,eigval] = KotteStabilityInfo(neweqpt',model,pvec);
    alleigval(:,ipt) = eigval;
    
    % set initial value loop (from different solution branches ?)
    iguess = alleqsol(:,ipt);
end
contpt = alleqsol; % (:,nvar*npar-(nvar-1):nvar*npar);   
varargout{1} = alleqflux;
varargout{2} = alleigval; % (:,nvar*npar-(nvar-1):nvar*npar);
% varargout{2} = eqpttype(:,end);
% all value
% varargout{3} = alleigval;
% varargout{4} = eqpttype;