function [contpt,varargout] = continuation(model,pvec,lambda,ap,iguess,opts)

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
nvar = size(iguess,2);
npts = size(iguess,1);
npar = length(lambda);
alleqsol = zeros(npts,npar*nvar);
eqpttype = zeros(npts,npar);
alleigval = zeros(nvar,npar*npts);
% continuation algorithm
% set parameter value loop
for ip = 1:npar 
    model.PM(ac-nvar) = lambda(ip);
    pvec(ap) = lambda(ip);
    sysfun = @(x)Kotte_givenNLAE(x,model,pvec);    
    % re-initialize neweqpt
    neweqpt= zeros(npts,nvar);
    for jval = 1:npts
        neweqpt(jval,:) = fsolve(sysfun,iguess(jval,:)',opts);
    end
    % store neweqpt
    alleqsol(:,(ip-1)*nvar+1:ip*nvar) = neweqpt;
    
    % get stability info
    [type,eigval] = KotteStabilityInfo(neweqpt,model,pvec);
    eqpttype(:,ip) = type;
    alleigval(:,(ip-1)*npts+1:ip*npts) = eigval;
    
    % set initial value loop (from different solution branches ?)
    iguess = neweqpt;
end
% last points only
contpt = alleqsol(:,nvar*npar-(nvar-1):nvar*npar);   
varargout{1} = alleigval(:,nvar*npar-(nvar-1):nvar*npar);
varargout{2} = eqpttype(:,end);
% all value
varargout{3} = alleigval;
varargout{4} = eqpttype;