% objective to minimize the L2-norm of fluxes (quadratic?)
function fx = objnoisy(x,p,data)

% assign parameter values from variable vector
% p(data.p_id) = x(data.nc+1:end);

% calculate flux
% vest = getmodelflux(x(1:data.nc),p);
% vest(1) = x(5).*x(3).*p(17)./(p(17)+x(4));

% calc objective
fx = sqrt(sum((x(data.nvar-data.nf*data.npert+1:data.nvar)-data.vexp).^2));