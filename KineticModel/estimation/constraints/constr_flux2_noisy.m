% objective function for estimation of kientic paramaters from noisy
% perturbation data
% concentrations are also variables
function fx_nlcon =...
        constr_flux2_noisy(x,p,p_id,ss_val)

ncons = 2*3+1;
np = size(ss_val,2);
% np - number of perturbation whose data is used
% fx_nlcon = zeros(ncons,1);

% parameters
K1ac = p(1);    % or 0.02
K3fdp = p(2);
L3fbp = p(3);
K3pep = p(4);
K2pep = p(5);
vemax = p(6);        % for bifurcation analysis: 0.7:0.1:1.3
KeFDP = p(7);        % or 0.45
ne = p(8);             % or 2
d = p(9);
V4max = p(10);
k1cat = p(11);   
V3max = p(12);    
V2max = p(13);  
ac = p(17);

% p(p_id) = x(4:5);
xss = ss_val(1:3,:);
fss = ss_val(end,:);
% model_flux = zeros(6,1);

% fluxes
model_flux(1) = k1cat.*x(3).*ac./(ac+K1ac) + x(end);
model_flux(2) = vemax.*(1-1./(1+(KeFDP./x(2)).^ne));
ratio = 1+x(2)./K3fdp;
model_flux(3) = V3max.*(ratio-1).*(ratio).^3./...
            (ratio.^4+L3fbp.*(1+x(1)./K3pep).^(-4));
model_flux(4) = x(5).*x(1)./(x(1)+x(4));    
model_flux(5) = V4max.*x(1);
model_flux(6) = d.*x(3);

% nl constraints
% constraint 1 : norm of flux for which parameters are estimated (1x1)
% model_flux = kotte_flux_noCAS(x(1:3),p);
% model_flux(1) = model_flux(1) + x(end);    
fx_nlcon(1,1) = sqrt(sum((model_flux(4)-fss).^2));

% constraint 2 : norm of all ss concentrations (mx1)
fx_nlcon(2:4,1) = sqrt(sum((repmat(x(1:3,1),1,np)-xss).^2,2));

% constraint 3 : ss for concentrations (mx1)
fx_nlcon(5:7,1) = abs(kotte_ode_eq(x(1:3),p,model_flux));

