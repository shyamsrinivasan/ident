function contpt = continuation(model,pvec,ap,iguess,opts)

global c

npts = size(iguess,1);
contpt = zeros(npts,size(iguess,2));

for i = 1:npts
    if sum(isnan(iguess(i,:)))==0
        pvec(ap) = 
        sysfun = @(x)Kotte_givenNLAE(x,model,pvec);
        [fpsol,fval,flag] = fsolve(sysfun,iguess(i,:),opts);
        fval = max(fval);
    else
        fpsol = NaN;
        fval = 1;
    end
    contpt(i,:) = fpsol;
    if exitflag~=1 || abs(fval)>opts
        contpt(i,:) = NaN;
    end
end

% continuation algorithm
% set parameter value loop
for ip = 1:npar
    lambda = alllambda(ip);
    model.PM(ap) = lambdal;
    pvec(ap) = lambda;
    sysfun = @(x)Kotte_givenNLAE(x,model,pvec);
    % set initial value loop (from different solution branches ?)
    iguess = neweqpt;
    % re-initialize neweqpt
    for jval = 1:nval
        neweqpt(jval,:) = fsolve(sysfun,iguess(jval,:),opts);
    end
    % store neweqpt
end
    