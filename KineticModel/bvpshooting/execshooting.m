function [yi,yf] = execshooting(fh,yinit,yterm,ti,tf,yiunkwn,yfknwn,delyi,delyf,yf,eps)
% execute shooting iteration
if nargin<11
    eps = 1e-4;
end

opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
[yi,yf,delyf] = itershooting(fh,yinit,yterm,ti,tf,yiunkwn,yfknwn,delyi,delyf,yf,opts);

while abs(delyf(yfknwn))>repmat(eps,length(yfknwn),1)
    [yi,yf,delyf] = itershooting(fh,yi,yf,ti,tf,yiunkwn,yfknwn,delyi,delyf,yf,opts);
end


