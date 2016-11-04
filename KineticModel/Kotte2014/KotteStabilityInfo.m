function [type,lambda] = KotteStabilityInfo(eqpts,model,pvec)
% type
% 1 : saddle 2D stable manifold
% 2 : saddle 1D stable manifold
% 3 : sink 3D stable manifold
% 4 : source 3D unstable manifold

% eigen values
lambda = zeros(size(eqpts,2),size(eqpts,1));
for ipts = 1:size(eqpts,1)
    [~,lambda(:,ipts)] = getKotteJacobian(eqpts(ipts,:)',pvec,model);
end
[neig,nvec] = size(lambda);

type = zeros(nvec,1);

eigval = real(lambda);

eigval(eigval>0) = 0;
eigval(eigval<0) = 1;

sumall = sum(eigval);

unstable = find(~all(eigval));
nunstable = length(unstable);

stable = find(all(eigval));
nstable = length(stable);

nnegeig_unstable = sumall(unstable);
nnegeig_stable = sumall(stable);

type(unstable(nnegeig_unstable==2)) = 1;
type(unstable(nnegeig_unstable==1)) = 2; 
type(unstable(nnegeig_unstable==0)) = 4; 

type(stable(nnegeig_stable==neig)) = 3;






