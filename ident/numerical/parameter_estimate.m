function flux3 = parameter_estimate(y, p)
	% function to solve for parameters of flux v3
	% input is experimental data
	% y - initial parameter estimates
	% p - experimental data [x1, x2, x3, ac, v1, v2, v3, v4]
	% convenience kinetics for flux 3
    fdp_sat = y(2) / K3fdp;
    pep_sat = y(1) / K3pep;
    nr_3 = V3max*fdp_sat;
    dr_3 = 1 + fdp_sat;
    regulation_activate = 1/(1 + 1/pep_sat);
    % regulation_inhibition = 1/(1 + pep_sat) for future reference
    flux3 = regulation_activate * nr_3/dr_3;