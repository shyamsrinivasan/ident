function [saveinterpt,status] =...
         getSintsecpt(USpoint,Seigw,model,pvec,opts,tfi,eps,alpha,itermax)
if nargin<9
    itermax = 100;
end
if nargin<8
    alpha = 0.8;
end
if nargin<7
    eps = 1e-3;
end
if nargin<6
    tfi = [1,5];
end

ip_interS = [USpoint+eps*Seigw USpoint-eps*Seigw];
saveeps = zeros(size(ip_interS,2),1);
saveinterpt = zeros(size(USpoint,1),2);
status = ones(size(ip_interS,2),1);

for ip = 1:size(ip_interS,2)
    xfwi = solveODEonly(1,ip_interS(:,ip),model,pvec,opts,tfi);
    iter = 1;
    while any(any(abs(xfwi-repmat(USpoint,1,size(xfwi,2)))>eps)) && iter<=itermax
        eps = eps*alpha;
        if ip==1
            ip_interS(:,ip) = USpoint+eps*Seigw;
        elseif ip==2
            ip_interS(:,ip) = USpoint-eps*Seigw;
        end 
        xfwi = solveODEonly(1,ip_interS(:,ip),model,pvec,opts,tfi);
        iter = iter+1;
    end
    if iter>itermax
        status(ip) = 0;
    end
    saveeps(ip) = eps;
    saveinterpt(:,ip) = ip_interS(:,ip);
end