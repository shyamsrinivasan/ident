function allsmplepar = sampleEKP(pvec,lb,ub,ap,npts)

nps = length(ap);
Nxeq = zeros(nps,npts);
for ip = 1:nps
    Nxeq(ip,:) = ivalueperturbation(lb(ap(ip)),ub(ap(ip)),npts);
end

allsmplepar = repmat(pvec,npts,1);
allsmplepar(:,ap) = Nxeq(:,:)';


