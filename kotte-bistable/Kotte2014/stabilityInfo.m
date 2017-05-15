function [eigvals_all,eigw] = stabilityInfo(funhandle,eqpts,model,pvec)

% eigen values
nvar = size(eqpts,2);
npts = size(eqpts,1);
eigvals_all = zeros(nvar,npts);
eigw = zeros(nvar*npts,nvar);
for ipts = 1:size(eqpts,1)
    [~,eigvals_all(:,ipts),w] =...
    getKotteJacobian(funhandle,eqpts(ipts,:)',pvec,model);
    % convert eigen vectors to unit eigen vectors  
    eigw((ipts-1)*nvar+1:ipts*nvar,:) = w./repmat(sqrt(sum(w.^2)),nvar,1);
end