function flux = Kotte_glycolysis_NLAE_flux(M,pvec,flux)
if nargin < 3
    flux = zeros(4,1);
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
    %parameters
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
    acetate = pvec(12);        % a.u acetate
end

%metabolic fluxes
%J(E, acetate)
flux(1) = kEcat.*M(1).*acetate./(acetate+KEacetate);

%vFbp(PEP,FBP)
ratio = 1+M(3)/KFbpFBP;
flux(3) = vFbpmax.*(ratio-1).*(ratio).^4/(ratio.^4+Lfbp*(1+M(2)./KFbpPEP));

%vEX(PEP)
flux(2) = vEXmax.*M(2)./(M(2)+KEXPEP);

%enzyme production fluxes
%E(FBP) for J (%FBP ---| Cra and Cra ---> E)
flux(4) = vemax.*(1-1./(1+(KeFBP./M(3)).^ne));
    