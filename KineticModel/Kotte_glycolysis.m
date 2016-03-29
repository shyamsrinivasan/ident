%model for bistability in gluconeogenesis, Kotte 2014
function dM = Kotte_glycolysis(t,M)
dM = zeros(4,1);

%parameters
kEcat
KEacetate
KFbpFBP
VFbpmax
Lfbp
KFbpPEP
vEXmax
KEXPEP
vemax
KeFBP
ne
d

%acetate --E--> PEP --vEX--> FBP --Fbp--> Nothing
% u(1) --M(1)--> M(2) --vEX--> M(3) --Fbp--> Nothing
%regulation
%PEP ---> Fbp
%FBP ---| Cra
%Cra ---> E

flux = Kotte_glycolysisflux(M);

%differential equations
%enzymes
%E
dM(1) = flux(4) - d*M(1);

%PEP
dM(2) = flux(1) - flux(2);

%FBP
dM(3) = flux(2) - flux(3);

%acetate
dM(4) = 0;

