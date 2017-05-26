function flux = kotte_convkinflux_noCAS(x,p)

npts = size(x,2);
flux = zeros(5,npts);

% parameters
K1ac = p(1);    % or 0.02
K3fdp = p(2);
% L3fdp = p(3);
K3pep = p(4);
K2pep = p(5);
vemax = p(6);        % for bifurcation analysis: 0.7:0.1:1.3
KeFDP = p(7);        % or 0.45
ne = p(8);             % or 2
% ac = p(9);
V4max = p(10);
k1cat = p(11);   
V3max = p(12);    
V2max = p(13);
K1pep = p(14);  % .02
K2fdp = p(15); % .3
rhoA = p(16); % 0.5

% convinience kientics for acetate uptake and conversion to pep
flux(1,:) = k1cat.*x(3,:).*(x(4,:)./K1ac)/...
            (1 + x(4,:)./K1ac + x(1,:)./K1pep);

% hill kinetics for enzyme production        
flux(2,:) = vemax.*(1-1./(1+(KeFDP./x(2,:)).^ne));

% convinience kinetics w/ allostery for vFbp 
acratio = x(1)./K3pep;
acflx = rhoA + (1-rhoA).*(acratio./(1+acratio));
% acflx = (rhoA + (1-rhoA).*acratio./(1+acratio)).^4;
flux(3,:) = V3max.*acflx.*(x(2,:)./K3fdp)./(1+x(2,:)./K3fdp);

% convinience kinetics for vEX(pep)
flux(4,:) = V2max.*(x(1,:)./K2pep)/(1+x(1,:)./K2pep+x(2,:)./K2fdp);

% linear kinetics for vPEPout
flux(5,:) = V4max.*x(1,:);