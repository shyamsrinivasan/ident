function flux = Kotte_glycolysisflux(M,flux)
if nargin < 2
    flux = zeros(4,1);
end

%acetate --E--> PEP --vEX--> FBP --Fbp--> Nothing
% u(1) or M(5) --M(1),flux(1)--> M(2) --vEX,flux(2)--> M(3) --Fbp,flux(3)--> Nothing
% E flux = flux(4)
%regulation
%PEP ---> Fbp
%FBP ---| Cra
%Cra ---> E

%metabolic fluxes
%J(E, acetate)
flux(1) = kEcat*M(1)*M(4)/(M(4)+KEacetate);

%vFbp(PEP,FBP)
ratio = 1+M(3)/KFbpFBP;
flux(3) = vFbpmax*(ratio-1)*(ratio)^4/(ratio^4+Lfbp*(1+M(2)/KFbpPEP));

%vEX(PEP)
flux(2) = vEXmax*M(2)/(M(2)+KEXPEP);

%enzyme production fluxes
%E(FBP) for J (%FBP ---| Cra and Cra ---> E)
flux(4) = vemax*(1-1/(1+(KeFBP/M(3))^ne));