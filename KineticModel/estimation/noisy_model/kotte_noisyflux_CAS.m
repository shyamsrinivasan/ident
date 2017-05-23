function [flux,FXflx] = kotte_noisyflux_CAS(x,p)

flux = cell(5,1);

% parameters
KEacetate = p(1);    % or 0.02
KFbpFBP = p(2);
Lfbp = p(3);
KFbpPEP = p(4);
KEXPEP = p(5);
vemax = p(6);        % for bifurcation analysis: 0.7:0.1:1.3
KeFBP = p(7);        % or 0.45
ne = p(8);             % or 2
acetate = p(9);
kPEPout = p(11);
kEcat = p(12);   
vFbpmax = p(13);    
vEXmax = p(14);   

noise = casadi.SX(rand(5,1)*2);

% metabolic fluxes
% J(E, acetate)
flux{1} = kEcat.*x(3).*acetate./(acetate+KEacetate) + noise(1);

% enzyme production fluxes
% E(FBP) for J (%FBP ---| Cra and Cra ---> E)
flux{2} = vemax.*(1-1./(1+(KeFBP./x(2)).^ne)) +noise(2);

% vFbp(PEP,FBP)
ratio = 1+x(2)./KFbpFBP;
flux{3} = vFbpmax.*(ratio-1).*(ratio).^3./...
            (ratio.^4+Lfbp.*(1+x(1)./KFbpPEP).^(-4)) + noise(3);

% vEX(PEP)
flux{4} = vEXmax.*x(1)./(x(1)+KEXPEP) + noise(4);

% vPEPout
flux{5} = kPEPout*x(1) + noise(5);

FXflx = casadi.Function('FXflx',{x,p},{[flux{1};...
                                        flux{2};...
                                        flux{3};...
                                        flux{4};...
                                        flux{5}]});