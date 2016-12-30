function [xdynr1,xeq] = calc1DWus(saddle,w,lambda,model,pvec,opts,tspanr,tspanf,eps)
% calculate 1D unstable manifold

w = w(:,real(lambda)>=0);

ival = [saddle+eps*w saddle-eps*w];
% choose either the positive or the negative perturbation for the manifold

xdynr1 = solveODEonly(1,ival(:,1),model,pvec,opts,tspanr);

% get eq points through perturbation and forward integration
[~,xeq1] = solveODEonly(1,ival(:,1),model,pvec,opts,tspanf);
[~,xeq2] = solveODEonly(1,ival(:,2),model,pvec,opts,tspanf);

xeq = [xeq1 xeq2];


