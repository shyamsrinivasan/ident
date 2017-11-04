% new contraints including eps in bounds
function fx = cons_typec_f3(x,p,data)

nvar = data.nvar;
nf = data.nf;
nc = data.nc;
npert = data.npert;
np = data.np;

pep = x(1:3:nc*npert);
fdp = x(2:3:nc*npert);
% E = x(3:3:nc*npert);
par = x(nc*npert+1:nc*npert+np);
flux = x(nc*npert+np+1:nc*npert+np+nf*npert);

ratio = 1+fdp./par(1);
fx1 = par(3).*(ratio-1).*(ratio).^3 - flux.*(ratio.^4+p(3).*(1+pep./par(2)).^(-4));

% concentration noise cons
fx2 = x(1:nc*npert) - data.xexp*(1+x(nvar-1));
fx3 = -x(1:nc*npert) + data.xexp*(1-x(nvar-1));

% flux noise cons
fx4 = flux - data.vexp*(1+x(nvar));
fx5 = -flux + data.vexp*(1-x(nvar));

fx = [fx1;fx2;fx3;fx4;fx5];
