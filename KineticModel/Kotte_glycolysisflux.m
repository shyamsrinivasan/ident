function flux = Kotte_glycolysisflux(M,flux)
if nargin < 2
    flux = zeros(4,1);
end
%parameters
kEcat = 1;
KEacetate = 0.1;    % or 0.02
KFbpFBP = 0.1;
vFbpmax = 1;
Lfbp = 4e6;
KFbpPEP = 0.1;
vEXmax = 1;
KEXPEP = 0.3;
vemax = 1.1;        % for bifurcation analysis: 0.7:0.1:1.3
KeFBP = 0.1;        % or 0.45
ne = 1;             % or 2

%acetate --E--> PEP --vEX--> FBP --Fbp--> Nothing
% u(1) or M(4) --M(1),flux(1)--> M(2) --vEX,flux(2)--> M(3) --Fbp,flux(3)--> Nothing
% E flux = flux(4)
%regulation
%PEP ---> Fbp
%FBP ---| Cra
%Cra ---> E

%metabolic fluxes
%J(E, acetate)
flux(1) = kEcat.*M(1).*M(4)./(M(4)+KEacetate);

%vFbp(PEP,FBP)
ratio = 1+M(3)/KFbpFBP;
flux(3) = vFbpmax.*(ratio-1).*(ratio).^4/(ratio.^4+Lfbp*(1+M(2)./KFbpPEP));

%vEX(PEP)
flux(2) = vEXmax.*M(2)./(M(2)+KEXPEP);

%enzyme production fluxes
%E(FBP) for J (%FBP ---| Cra and Cra ---> E)
flux(4) = vemax.*(1-1./(1+(KeFBP./M(3)).^ne));