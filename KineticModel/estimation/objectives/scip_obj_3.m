function fx = scip_obj_3(x,p)

% x = [K3pep;V3max;K3fdp;L3fdp];
% p = [pep;fdp;enz;ac;gflux];

ratio = 1+p(2)./x(3);
flux = x(2).*(ratio-1).*(ratio).^3./...
            (ratio.^4+x(4).*(1+p(1)./x(1)).^(-4));
fx = sqrt((flux-p(5))^2);