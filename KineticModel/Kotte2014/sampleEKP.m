function allsmplepar = sampleEKP(pvec,lb,ub,idp,npts)

nps = length(idp);
Nxeq = zeros(nps,npts);
for ip = 1:nps
    if lb(ip)<ub(ip)
        Nxeq(ip,:) = ivalueperturbation(lb(ip),ub(ip),npts);
    else
        % ub becomes lb and lb becomes ub
        Nxeq(ip,:) = ivalueperturbation(ub(ip),lb(ip),npts);
    end
end

allsmplepar = repmat(pvec,npts,1);
allsmplepar(:,idp) = Nxeq(:,:)';


