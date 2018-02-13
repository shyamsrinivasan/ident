function dy = kotteCKNLAE(y, p)
	dy = zeros(3, 1);
	flux = kotteCKflux(y, p);

	dy(1) = flux(1) - flux(4) - flux(5);
    dy(2) = flux(4) - flux(3);
    dy(3) = flux(2) - flux(6);
