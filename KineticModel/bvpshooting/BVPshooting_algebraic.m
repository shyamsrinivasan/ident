function fx = BVPshooting_algebraic(delyi,delyf,xic,yiknwn,yfknwn)
% x - n x n-r matrix of solutions from n-r reverse integration runs on the
%             adjoint equation
% 
% generate the n-r algebraic equations to be solved to get n-r unknown
% initial conditions
% xr+1*delyr+1 + ... + xn*delyn
r = length(yiknwn);
nvar = size(delyi,1);
yiunkwn = setdiff(1:nvar,yiknwn);

fx = zeros(nvar-r,1);
for m = 1:nvar-r
    fx(m) = sum(xic(yiunkwn,m).*delyi(yiunkwn))-delyf(yfknwn(m));
end


% for m = 1:nvar-r
%     fx(m) = sum(xic(yfknwn,m).*delyi(yfknwn))-delyf(yfknwn(m));
% end
