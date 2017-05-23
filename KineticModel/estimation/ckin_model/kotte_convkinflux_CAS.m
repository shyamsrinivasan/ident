function [flux,FXflx] = kotte_convkinflux_CAS(x,p)

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
K1pep = p(15);  % .02
KEXfdp = p(16); % .3
rhoA = p(17); % 0.5

% convinience kientics for acetate uptake and conversion to pep
flux{1} = kEcat.*x(3).*(acetate./KEacetate)/...
            (1 + acetate./KEacetate + x(1)./K1pep);

% hill kinetics for enzyme production        
flux{2} = vemax.*(1-1./(1+(KeFBP./x(2)).^ne));

% convinience kinetics w/ allostery for vFbp 
acratio = x(1)./KFbpPEP;
acflx = (rhoA + (1-rhoA).*acratio./(1+acratio)).^4;
flux{3} = vFbpmax.*acflx.*(x(2)./KFbpFBP)./(1+x(2)./KFbpFBP);

% convinience kinetics for vEX(pep)
flux{4} = vEXmax.*(x(1)./KEXPEP)/(1+x(1)./KEXPEP+x(2)./KEXfdp);

% linear kinetics for vPEPout
flux{5} = kPEPout*x(1);

FXflx = casadi.Function('FXflx',{x,p},{[flux{1};...
                                        flux{2};...
                                        flux{3};...
                                        flux{4};...
                                        flux{5}]});

