function constr = contr_flux1_noisy(x,p,p_id,gflux,gcon)

% x = [pep;fdp;enz;ac;K1ac;k1cat];

% calculate flux
p(p_id) = x(5:end);
flux = kotte_flux(x(1:4),p);

% constr : noisy flux error
noisy_flux = flux + e;

% constr 1 : flux norm for optimized parameters
constr1 = norm(noisy_flux(1)-gflux);

% constr 2 : concentration steady state
constr2 = kotte_ode(x(1:4),p,noisy_flux);

% constr 3 : concentration norm
constr3 = norm(x(1:4)-gcon);

constr = [constr1;constr2;constr3];

function dx = kotte_ode(x,p,flux)
dx = zeros(4,1);

dx(1) = flux(1)-flux(4)-flux(5);
dx(2) = flux(4)-flux(3);
dx(3) = flux(2)-p(9)*x(3);
dx(4) = x(4)-x(4);



