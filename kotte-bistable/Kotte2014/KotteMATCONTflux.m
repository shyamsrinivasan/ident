function flux = KotteMATCONTflux(M,pvec,flux)
if nargin < 3
    flux = zeros(5,1);
end

if nargin < 2    
    fprintf('No parameter vector\n');
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
    acetate = 2;        % a.u acetate
    kPEPout = 0.2;
else
    % parameters    
    KEacetate = pvec(1);    % or 0.02
    KFbpFBP = pvec(2);
    Lfbp = pvec(3);
    KFbpPEP = pvec(4);
    KEXPEP = pvec(5);
    vemax = pvec(6);        % for bifurcation analysis: 0.7:0.1:1.3
    KeFBP = pvec(7);        % or 0.45
    ne = pvec(8);             % or 2
    acetate = pvec(9);
    kPEPout = pvec(11);
    kEcat = pvec(12);   
    vFbpmax = pvec(13);    
    vEXmax = pvec(14);   
end

%acetate --E--> PEP --vEX--> FBP --Fbp--> Nothing
% u(1) or M(4) --M(1),flux(1)--> M(2) --vEX,flux(2)--> M(3) --Fbp,flux(3)--> Nothing
% E flux = flux(4)
%regulation
%PEP ---> Fbp
%FBP ---| Cra
%Cra ---> E

% metabolic fluxes
% J(E, acetate)
flux(1) = kEcat.*M(3).*acetate./(acetate+KEacetate);

% enzyme production fluxes
% E(FBP) for J (%FBP ---| Cra and Cra ---> E)
flux(2) = vemax.*(1-1./(1+(KeFBP./M(2)).^ne));

% vFbp(PEP,FBP)
ratio = 1+M(2)/KFbpFBP;
flux(3) = vFbpmax.*(ratio-1).*(ratio).^3/(ratio.^4+Lfbp*(1+M(1)./KFbpPEP).^(-4));

% vEX(PEP)
flux(4) = vEXmax.*M(1)./(M(1)+KEXPEP);

% vPEPout(PEP)
flux(5) = kPEPout*M(1);



