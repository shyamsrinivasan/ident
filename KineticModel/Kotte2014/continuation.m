function contpt = continuation(systemfun,iguess,opts)

global c

npts = size(iguess,1);
contpt = zeros(npts,size(iguess,2));

for i = 1:npts
    if sum(isnan(iguess(i,:)))==0
        [fpsol,fval,flag] = fsolve(systemfun,iguess(i,:),opts);
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