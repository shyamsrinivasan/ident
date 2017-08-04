% objective to minimize weighted sum of L2-norm of fluxes and
% concentrations
function fx = obj_typeb(x,p,data)

nc = data.nc;
nf = data.nf;
npert = data.npert;
nvar = data.nvar;
w1 = data.flux_wt;
w2 = data.conc_wt;

% flux norm
flux_norm = sqrt(sum((x(nvar-nf*npert+1:nvar)-data.vexp).^2));

% concentration norm
conc_norm = sqrt(sum((x(1:nc*npert)-data.xexp).^2));

% calc objective
fx = w1.*flux_norm + w2.*conc_norm;

