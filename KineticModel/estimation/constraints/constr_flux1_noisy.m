function nlcon = constr_flux1_noisy(x,optim_p,p_id)

% optimp_p = [flux parameters(w/o estimated parameters);xss;fss]

% all options passed as options structure
% if isfield(optim_p,'p')
%     optim_p = optim_p.p;
% end
% if isfield(optim_p,'p_id')
%     p_id = optim_p.p_id;
% end
% if isfield(optim_p,'xss')
%     xss = optim_p.xss;
% end
% if isfield(optim_p,'fss')
%     fss = optim_p.fss;
% end
fss = optim_p(end);
xss = optim_p(end-4:end-1);
optim_p(end-4:end) = [];

% variables
% x = [pep;fdp;enz;ac;K1ac;k1cat;e];
optim_p(p_id) = x(5:4+length(p_id));

% nl constraints
% constraint 1 : norm of flux for which parameters are estimated (1x1)
model_flux = kotte_flux_noCAS(x(1:4),optim_p);
model_flux(1) = model_flux(1) + x(end);    
nlcon1 = sqrt(sum((model_flux(1)-fss).^2));

% constraint 2 : norm of all ss concentrations (1x1)
nlcon2 = sqrt(sum((x(1:4)-xss).^2));

% constraint 3 : ss for concentrations (mx1)
nlcon3 = kotte_ode(x(1:4),optim_p,model_flux);

% all constraints : (m+2)x1
nlcon = [nlcon1;nlcon2;nlcon3];

function dx = kotte_ode(x,p,flux)
dx = zeros(4,1);

dx(1) = flux(1)-flux(4)-flux(5);
dx(2) = flux(4)-flux(3);
dx(3) = flux(2)-p(9)*x(3);
dx(4) = x(4)-x(4);



