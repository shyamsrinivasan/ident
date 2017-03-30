% Figure 3 - separatrix based on stable manifold theorem and numerical
% approximation from Seydel 2009 and Lyons et al., 2014
runKotte

% get saddle node
[saddle,saddlepar] = getsaddlenode(data.s1,data.x1,5e-3);
pvec(ap) = saddlepar;
model.PM(ac-length(saddle)) = saddlepar;

% perturb saddle to get steady states
eps = 1e-4;
tspanf = 0:0.1:2000;
pival = saddle+eps*[1;1;1];
[~,xeq1,~,feq1] = solveODEonly(1,pival,model,pvec,opts,tspanf);
nival = saddle-eps*[1;1;1];
[~,xeq2,~,feq2] = solveODEonly(1,nival,model,pvec,opts,tspanf);

% calculate and plot separatrix as 2-D projections
% stable 3D manifolds
tspanr = [0 -8.5];
NumericalSeparatrix(model,pvec,opts,ap,data.s1,data.x1,[xeq1 xeq2],'stable',tspanr,3,5e-3);

% stable 2D manifolds
tspanr = [0 -8.5];
NumericalSeparatrix(model,pvec,opts,ap,data.s1,data.x1,[xeq1 xeq2],'stable',tspanr,2,5e-3);

% unstable 3D manifolds
tspanr = [0 -20];
NumericalSeparatrix(model,pvec,opts,ap,data.s1,data.x1,[xeq1 xeq2],'unstable',tspanr,3,5e-3);

% unstable 2D manifolds
tspanr = [0 -20];
NumericalSeparatrix(model,pvec,opts,ap,data.s1,data.x1,[xeq1 xeq2],'unstable',tspanr,2,5e-3);