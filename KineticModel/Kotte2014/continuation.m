function contpt = continuation(model,pvec,lambda,ap,iguess,opts)

global c

npts = size(iguess,1);
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
    % alleqsol
    % set initial value loop (from different solution branches ?)
    iguess = neweqpt;
end
    