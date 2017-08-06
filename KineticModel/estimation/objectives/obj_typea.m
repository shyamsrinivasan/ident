% objective to minimize the L2-norm of fluxes (quadratic?)
function fx = obj_typea(x,p,data)

% calc objective
fx = sqrt(sum((x(data.nvar-data.nf*data.npert+1:data.nvar)-data.vexp).^2));