function [flux,FXflx] = kotte_flux_CAS(x,p)

flux = cell(6,1);

% parameters
K1ac = p(1);    % or 0.02
K3fdp = p(2);
L3fdp = p(3);
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

% noise = casadi.SX(rand(5,1)*2);

% metabolic fluxes
% J(E, acetate)
flux{1} = k1cat.*x(3).*ac./(ac+K1ac);

% enzyme production fluxes
% E(FBP) for J (%FBP ---| Cra and Cra ---> E)
flux{2} = vemax.*(1-1./(1+(KeFDP./x(2)).^ne));

% vFbp(PEP,FBP)
ratio = 1+x(2)./K3fdp;
flux{3} = V3max.*(ratio-1).*(ratio).^3./...
            (ratio.^4+L3fdp.*(1+x(1)./K3pep).^(-4));

% vEX(PEP)
flux{4} = V2max.*x(1)./(x(1)+K2pep);

% vPEPout
flux{5} = V4max.*x(1);

% E --->
flux{6} = d.*x(3);

FXflx = casadi.Function('FXflx',{x,p},{[flux{1};...
                                        flux{2};...
                                        flux{3};...
                                        flux{4};...
                                        flux{5};...
                                        flux{6}]});