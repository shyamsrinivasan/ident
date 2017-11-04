% function to evaluate jacobians/gradients based on the given function
% using ADMAT
% Note: 
% Function and every sub-function within must use the 'cons' function
% from ADMAT to provide the correct jacobian.
% A forward difference approximation is also availabel to be used without 
% using ADMAT and its 'cons' method
function Gradient(fhandle,xeq,method)
% fhandle - function handle to calculate jacobian/gradient
% xeq - point at which grad is to calculated
% method - finite difference or {ADMAT (AD)}

if strcmpi(method,'AD')
elseif strcmpi(method,'FD')
end

function AD()
return

function FD()
end