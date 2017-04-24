function [y] = bifur_direct()
% solve equations in bifurequations.m using Newton's method
% requires 
% 1. jacobian of all equations in bifurequations.m

% inputs
% lambdai - parameter vector at each iteration

% y = [state variables,x;...
%      power margin,l;...
%      left eigen vector,w];

% solve the following set of equations
% f(x,lambda1) = 0
% w1*Dfx = 0
% w1c-1 = 0


lambdai