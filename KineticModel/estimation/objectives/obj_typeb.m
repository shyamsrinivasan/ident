% objective to minimize weighted sum of L2-norm of fluxes and
% concentrations
function fx = obj_typeb(x,p,data)

nc = data.nc;
nf = data.nf;
npert = data.npert;
nvar = data.nvar;
np = data.np;
w1 = data.flux_wt;
w2 = data.conc_wt;

% flux norm
flux_norm = sqrt(sum((x(nc*npert+np+1:nc*npert+np+nf*npert)-data.vexp).^2));

% concentration norm
conc_norm = sqrt(sum((x(1:nc*npert)-data.xexp).^2));

% calc objective
fx = w1.*flux_norm + w2.*conc_norm + x(nvar-1).^4 + x(nvar).^4;

