function [yi,yf,delyf] =...
execshooting(fh,yi,yterm,ti,tf,yiunkwn,yfknwn,delyi,delyf,yf,eps,opts,niter)
% execute shooting iteration
if nargin<13
    niter = 100;
end
if nargin<12
    opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
end
if nargin<11
    eps = 1e-4;
end
if nargin<10||isempty(yf)
    yf = zeros(size(yi,1),1);
end
if nargin<9||isempty(delyf)
    delyf = zeros(size(yi,1),1);
end

iter = 1;
while abs(delyf(yfknwn))>repmat(eps,length(yfknwn),1) 
    fprintf('Iteration #%d\n',iter);
    if iter<niter
        [yi,yf,delyf] = itershooting(fh,yi,yf,ti,tf,yiunkwn,yfknwn,delyi,delyf,yf,opts);
        iter = iter+1;
    end
end


