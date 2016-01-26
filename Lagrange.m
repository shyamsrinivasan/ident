function [L,Ldot] = Lagrange(tau,j,x)
%Lagrange polynomial - L(tau)
K = size(x,1)-1;
idx1 = setdiff(0:K, j);
nr_x1 = repmat(tau,length(idx1),1)-x(idx1+1);
dr_x1 = repmat(x(j+1),length(idx1),1) - x(idx1+1);
L = prod(nr_x1./dr_x1);

nr_t = 0;
dr_t = 1;
k = setdiff(0:K,j);
for jk = 1:length(k)
    idx = setdiff(0:K, [k(jk),j]);
    delx = repmat(tau,length(idx+1),1)-x(idx+1);
    nr_t = nr_t + prod(delx);
    dr_t = dr_t*(x(j+1)-x(k(jk)+1));
end
Ldot = nr_t/dr_t;
