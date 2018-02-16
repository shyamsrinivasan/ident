function residual = parameter_estimate(p, y)
	% function to solve for parameters of flux v3
	% input is experimental data
	% p - initial parameter estimates p = [V3max, K3fdp, K3pep]
	% y - experimental data [ac1, x11, x21, x31, v11, v21, v31, v41,...
    %                        ac2, x12, x22, x32, v12, v22, v32, v42,...
    %                        ac3, x13, x23, x33, v13, v23, v33, v43]
	% form nonlinear algebraic equations
    residual = zeros(3, 1);
    
    % get flux from input experimental data
    flux = y(7:8:end)';
    
    % calculate rhs of flux equation based on input parameters in p
    rhs = flux3_CK(y, p);
    
    residual(:) = flux - rhs';    
return

function flux = flux3_CK(y, p)
% calculate value of flux v3 - fbp in kotte model using MWC model

fdp_sat = y(3:8:end)./p(2);
pep_sat = y(1:8:end)./p(3);
nr_3 = p(1).*fdp_sat;
dr_3 = 1 + fdp_sat;
regulation_activate = 1./(1 + 1./pep_sat);
% regulation_inhibition = 1/(1 + pep_sat) for future reference
flux = regulation_activate.*nr_3./dr_3;

return