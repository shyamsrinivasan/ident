function fx = consnoisyf2(x,p,data)

nvar = data.nvar;
nf = data.nf;
nc = data.nc;
npert = data.npert;
p_id = data.p_id;

pep = x(1:3:nc*npert);
% fdp = x(2:3:nc*npert);
% E = x(3:3:nc*npert);
par = x(nc*npert+1:nc*npert+length(p_id));
flux = x(nvar-nf*npert+1:nvar);

fx = par(2).*pep - flux.*(pep+par(1));