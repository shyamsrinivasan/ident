function allsmplepar = sampleEKP(pvec,lb,ub,idp,npts)

nps = length(idp);
Nxeq = zeros(nps,npts);
for ip = 1:nps
    Nxeq(ip,:) = ivalueperturbation(lb(ip),ub(ip),npts);
end

allsmplepar = repmat(pvec,npts,1);
allsmplepar(:,idp) = Nxeq(:,:)';


