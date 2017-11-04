function fx = cons_typeb_f3(x,p,data)

nvar = data.nvar;
nf = data.nf;
nc = data.nc;
npert = data.npert;
p_id = data.p_id;

pep = x(1:3:nc*npert);
fdp = x(2:3:nc*npert);
% E = x(3:3:nc*npert);
par = x(nc*npert+1:nc*npert+length(p_id));
flux = x(nvar-nf*npert+1:nvar);

ratio = 1+fdp./par(1);
fx = par(3).*(ratio-1).*(ratio).^3 - flux.*(ratio.^4+p(3).*(1+pep./par(2)).^(-4));