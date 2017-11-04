% design algorithm for bifurcation studies
% build and test model to be used for design
% runKotte

% iterative method to calculate closest bifurcation
% 1. get initial direction of parameter search n0
% 2. compute saddle node bifurcations along n0
%   - solve the following system using Newton's or another method to get
%   li, lambdai, xi st lambdai = lambda0 + n0*li
%   f(x1,lambda0+ln0) = 0
%   w1*Df(x1,lambda0+ln0) = 0
%   w1c - 1 = 0
% where c is any constant vector of appropriate size to make sure w1 is
% nonzero
% w1 - left eigen vector
% 3. set ni = wi*flambda (sensitivty w.r.t parameters lambda)
% 4. iterate through steps 1,2 and 3.

% step 1. 
% get initial direction vector n0 - 0.2 for reactive power loads and 1.0
% for all real power loads.

% step 2. 
% solve system of nonlinear AE to get all necessary values

runKotte
[orig_saddle,orig_saddlepar] = getsaddlenode(data.s1,data.x1,5e-3);
pvec(ap) = orig_saddlepar;
model.PM(ac-length(orig_saddle)) = orig_saddlepar;

% perturb saddle to get steady states
eps = 1e-4;
tspanf = 0:0.1:2000;
pival = orig_saddle+eps*[1;1;1];
[~,xeq1,~,feq1] = solveODEonly(1,pival,model,pvec,opts,tspanf);
nival = orig_saddle-eps*[1;1;1];
[~,xeq2,~,feq2] = solveODEonly(1,nival,model,pvec,opts,tspanf);
xss = [xeq1 xeq2];

[J,eigval] = getKotteJacobian(@Kotte_givenNLAE,orig_saddle,pvec,model);
[v,d,w] = eig(J);

