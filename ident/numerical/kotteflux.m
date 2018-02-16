function flux = kotteflux(y, p)
	% all model parameters
	K1ac = p(1);
    K3fdp = p(2);
    L3fdp = p(3);
    K3pep = p(4);
    K2pep = p(5);
    vemax = p(6);
    Kefdp = p(7);
    ne = p(8);
    d = p(9);
    V4max = p(10);
    k1cat = p(11);
    V3max = p(12);
    V2max = p(13);
    ac = p(14);

    flux_1 = k1cat.*y(3, :).*ac./(ac + K1ac);
    flux_2 = vemax.*(1 - 1./(1 + (Kefdp./y(2, :)).^ne));
    % convenience kinetics for flux 3
    fdp_sat = 1 + y(2, :)./K3fdp;
    pep_sat = 1 + y(1, :)./K3pep;
    % v3 with MWC kinetics
    flux_3 = V3max.*(fdp_sat - 1).*(fdp_sat.^3)./...
            (fdp_sat.^4+L3fdp.*(pep_sat.^(-4)));
    flux_4 = V2max.*y(1, :)./(y(1, :) + K2pep);
    flux_5 = V4max.*y(1, :);
    flux_6 = d.*y(3, :);
    flux = [flux_1;flux_2;flux_3;flux_4;flux_5;flux_6]';

    
    