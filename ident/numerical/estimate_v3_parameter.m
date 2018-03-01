function residual = estimate_v3_parameter(p, y)
% function to solve for parameters of flux v3
	% input is experimental data
	% p - initial parameter estimates p = [V3max, K3fdp, K3pep]
	% y - experimental data [ac1, x11, x21, x31, v11, v21, v31, v41,...
    %                        ac2, x12, x22, x32, v12, v22, v32, v42,...
    %                        ac3, x13, x23, x33, v13, v23, v33, v43]
	% form nonlinear algebraic equations
    residual = zeros(3,1);    
    % calculate rhs of flux equation based on input parameters in p
    rhs = flux3_MWC(y, p);
    % get flux from input experimental data
    flux = y(7:8:end)';
    residual(:) = flux - rhs';     
return

function flux = flux3_MWC(y, p)
% calculate value of flux v3 - fbp in kotte model using MWC model

fdp_sat = 1 + y(3:8:end)./p(2);
pep_sat = 1 + y(1:8:end)./p(3);
flux = p(1).*(fdp_sat - 1).*(fdp_sat.^3)./...
       (fdp_sat.^4+p(4).*(pep_sat.^(-4)));
return