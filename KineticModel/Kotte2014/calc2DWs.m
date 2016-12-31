function [xchop,ychop,zchop] =...
         calc2DWs(saddle,saddlepar,ap,model,pvec,tspanr,opts)

ac = find(strcmpi(model.mets,'ac[e]'));     

% one single saddle node with 2D stable manifold (local Ws)
pvec(ap) = saddlepar;
model.PM(ac-length(saddle)) = saddlepar;

% compute Jacobian, eigenvalues, eigenvector at saddle points with 2D stable
% manifolds
[~,eig,eigvec] = getKotteJacobian(saddle,pvec,model);

% obtain the 2 eigenvectors of 2D linear eigenspace
stableeigvec = eigvec(:,eig<0);

% compute invariant manifold via numerical integration from a circular 
% locus in the stable eigen space from which reverse
% intergration is to be performed to identify the 2D stable manifold
% surface

% Lyons et al., 2014 code
% circle parameters
points = 801;
radius = 0.01;

% obtain coordinates of circle with radius r in (x1,x2) plane
[x1,x2] = getplanarcircle(points,radius);

% perform linear mapping of unit circle onto plane in R3 spanned by W1,W2
circlenew = manifoldlinearmapping(x1,x2,stableeigvec(:,1),stableeigvec(:,2));

% translate mapping to saddle point
circlenew = circlenew + repmat(saddle,1,size(circlenew,2));
circlenew = circlenew';

[x,y,z,dynr] = get2Dmanifoldpoints(circlenew,model,pvec,tspanr,opts);
% cut put values that are outside the area of interest
[xchop,ychop,zchop,r] = chopvals(x,y,z,[2.5 2.5 2.5]);

% [xnew,ynew,znew] = removeredundantpoints(real(xchop),real(ychop),real(zchop),0.01);
% Manifold2DPlot(real(xnew),real(ynew),real(znew));