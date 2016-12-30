function [xdynr1,xeq] = calc1DWs(saddle,w,lambda,model,pvec,opts,tspanr)
% calculate 1D stable manifold

w = w(:,real(lambda)<0);

ival = [saddle+eps*w saddle-eps*w];
% choose either the positive or the negative perturbation for the manifold
xdynr1 = solveODEonly(1,ival(:,1),model,pvec,opts,tspanr);

% get eq points through perturbation and forward integration
[~,xeq1] = solveODEonly(1,ival(:,1),model,pvec,opts,tspanf);
[~,xeq2] = solveODEonly(1,ival(:,2),model,pvec,opts,tspanf);

xeq = [xeq1 xeq2];