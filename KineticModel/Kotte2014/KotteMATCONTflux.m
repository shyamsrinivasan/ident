function flux = KotteMATCONTflux(M,pvec,flux)
if nargin < 4
    flux = zeros(4,1);
end
if nargin < 3
    model = struct([]);
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
else
    % parameters    
    kEcat = pvec(1);
    KEacetate = pvec(2);    % or 0.02
    KFbpFBP = pvec(3);
    vFbpmax = pvec(4);
    Lfbp = pvec(5);
    KFbpPEP = pvec(6);
    vEXmax = pvec(7);
    KEXPEP = pvec(8);
    vemax = pvec(9);        % for bifurcation analysis: 0.7:0.1:1.3
    KeFBP = pvec(10);        % or 0.45
    ne = pvec(11);             % or 2
    acetate = pvec(12);
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



