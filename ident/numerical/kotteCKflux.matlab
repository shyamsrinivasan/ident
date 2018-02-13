function flux = kotteCKflux(y, p)
	% all model parameters
	K1ac = p(1)
    K3fdp = p(2)
    K3pep = p(3)
    K2pep = p(4)
    vemax = p(5)
    Kefdp = p(6)
    ne = p(7)
    d = p(8)
    V4max = p(9)
    k1cat = p(10)
    V3max = p(11)
    V2max = p(12)
    ac = p(13)

    flux_1 = k1cat * y(3) * ac / (ac + K1ac)
    flux_2 = vemax * (1 - 1 / (1 + (Kefdp / y(2))^ne))
    % convenience kinetics for flux 3
    fdp_sat = y(2) / K3fdp
    pep_sat = y(1) / K3pep
    nr_3 = V3max*fdp_sat
    dr_3 = 1 + fdp_sat
    regulation_activate = 1/(1 + 1/pep_sat)
    % regulation_inhibition = 1/(1 + pep_sat) for future reference
    flux_3 = regulation_activate * nr_3/dr_3
    flux_4 = V2max * y(1) / (y(1) + K2pep)
    flux_5 = V4max * y(1)
    flux_6 = d * y(3)
    flux = [flux_1, flux_2, flux_3, flux_4, flux_5, flux_6]'
