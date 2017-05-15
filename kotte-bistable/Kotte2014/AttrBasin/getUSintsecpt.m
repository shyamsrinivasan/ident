function [saveinterpt,status] =...
         getUSintsecpt(USpoint,USeigw,model,pvec,opts,trend,eps,alpha,itermax)
if nargin<9
    itermax = 100;
end
if nargin<8
    alpha = 0.8;
end
if nargin<7
    eps = 1e-3;
end

ip_interUS = [USpoint+eps*USeigw USpoint-eps*USeigw];
saveeps = zeros(size(ip_interUS,2),1);
saveinterpt = zeros(size(USpoint,1),2);
status = ones(size(ip_interUS,2),1);

for ip = 1:size(ip_interUS,2)
    xrev = solveODEonly(1,ip_interUS(:,ip),model,pvec,opts,[0,-trend]);
    iter = 1;
    while any(any(abs(xrev-repmat(USpoint,1,size(xrev,2)))>eps)) && iter<=itermax
        eps = eps*alpha;
        if ip==1
            ip_interUS(:,ip) = USpoint+eps*USeigw;
        elseif ip==2
            ip_interUS(:,ip) = USpoint-eps*USeigw;
        end 
        xrev = solveODEonly(1,ip_interUS(:,ip),model,pvec,opts,[0,-trend]);
        iter = iter+1;
    end
    if iter>itermax
        status(ip) = 0;
    end
    saveeps(ip) = eps;
    saveinterpt(:,ip) = ip_interUS(:,ip);
end

