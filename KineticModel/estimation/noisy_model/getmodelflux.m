function flux = getmodelflux(x,p)

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

% fluxes
flux(1) = k1cat.*x(3).*ac./(ac+K1ac);
flux(2) = vemax.*(1-1./(1+(KeFDP./x(2)).^ne));
ratio = 1+x(2)./K3fdp;
flux(3) = V3max.*(ratio-1).*(ratio).^3./...
            (ratio.^4+L3fbp.*(1+x(1)./K3pep).^(-4));
flux(4) = V2max.*x(1)./(x(1)+K2pep);    
flux(5) = V4max.*x(1);
flux(6) = d.*x(3);