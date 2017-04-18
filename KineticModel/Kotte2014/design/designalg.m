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

% 
