%model for bistability in gluconeogenesis, Kotte 2014
function dM = Kotte_glycolysis(t,M,pvec)
dM = zeros(3,1);

% parameters
d = pvec(13); % = 0.18; %or 0.25 or 0.35

%acetate --E--> PEP --vEX--> FBP --Fbp--> Nothing
% u(1) --M(1)--> M(2) --vEX--> M(3) --Fbp--> Nothing
%regulation
%PEP ---> Fbp
%FBP ---| Cra
%Cra ---> E

flux = Kotte_glycolysisflux(M,pvec);

%differential equations
%enzymes
%E
dM(1) = flux(4) - d*M(1);

%PEP
dM(2) = flux(1) - flux(2);

%FBP
dM(3) = flux(2) - flux(3);



