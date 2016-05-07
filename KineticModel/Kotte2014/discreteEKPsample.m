function allsmplepar = discreteEKPsample(pvec,vals,idp,npts)
vals = ToColumnVector(vals);

rndselect = randi(length(vals),[length(idp) npts]);
nps = length(idp);
Nxeq = zeros(nps,npts);
for ip = 1:npts
    Nxeq(:,ip) = vals(rndselect(:,ip))';
end

allsmplepar = repmat(pvec,npts,1);
allsmplepar(:,idp) = Nxeq(:,:)';