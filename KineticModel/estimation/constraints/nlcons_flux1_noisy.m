% nl constraints for parameter estimation
function nlcons_flux1_noisy(x,p,data)
% x = [x,pi];

% assign parameter values from variable vector
p(data.p_id) = x(data.nc+1:end);

% calculate flux
vest = getmodelflux(x(1:data.nc),p);


