% function [type,eigvals_all] = KotteStabilityInfo(eqpts,model,pvec)
% get stability information for all points passed as input
% requires ADMAT to be installed for algorithmic differentiation (AD) to
% get jacobian and corresponding eigen values at each eqpts passed in
function [type,eigvals_all,eigw] = KotteStabilityInfo(eqpts,model,pvec)
% type
% 1 : saddle 2D stable manifold
% 2 : saddle 1D stable manifold
% 3 : sink 3D stable manifold
% 4 : source 3D unstable manifold

% eigen values
nvar = size(eqpts,2);
npts = size(eqpts,1);
eigvals_all = zeros(nvar,npts);
eigw = zeros(nvar*npts,nvar);
for ipts = 1:size(eqpts,1)
    [~,eigvals_all(:,ipts),w] = getKotteJacobian(eqpts(ipts,:)',pvec,model);
    eigw((ipts-1)*nvar+1:ipts*nvar,:) = w;
end
[neig,nvec] = size(eigvals_all);
type = zeros(nvec,1);

eigval = real(eigvals_all);

eigval(eigval>0) = 0;
eigval(eigval<0) = 1;

sumall = sum(eigval);

unstable = find(~all(eigval));
% nunstable = length(unstable);
stable = find(all(eigval));
% nstable = length(stable);
nnegeig_unstable = sumall(unstable);
nnegeig_stable = sumall(stable);

type(unstable(nnegeig_unstable==2)) = 1;
type(unstable(nnegeig_unstable==1)) = 2; 
type(unstable(nnegeig_unstable==0)) = 4;
type(stable(nnegeig_stable==neig)) = 3;






