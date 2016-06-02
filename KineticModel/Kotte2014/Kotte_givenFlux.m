function flux = Kotte_givenFlux(M,pvec,model,flux)
if nargin < 4
    flux = zeros(5,1);
    flux = cons(flux,M);
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
    kPEPout = pvec(11);
    kEcat = pvec(12);   
    vFbpmax = pvec(13);    
    vEXmax = pvec(14);   
end

pep = strcmpi(model.mets,'pep[c]');
fdp = strcmpi(model.mets,'fdp[c]');
enz = strcmpi(model.mets,'enz[c]');
ac = strcmpi(model.mets,'ac[e]');


%acetate --E--> PEP --vEX--> FBP --Fbp--> Nothing
% u(1) or M(4) --M(1),flux(1)--> M(2) --vEX,flux(2)--> M(3) --Fbp,flux(3)--> Nothing
% E flux = flux(4)
%regulation
%PEP ---> Fbp
%FBP ---| Cra
%Cra ---> E

% metabolic fluxes
% J(E, acetate)
flux(1) = kEcat.*M(enz).*M(ac)./(M(ac)+KEacetate);

% enzyme production fluxes
% E(FBP) for J (%FBP ---| Cra and Cra ---> E)
flux(2) = vemax.*(1-1./(1+(KeFBP./M(fdp)).^ne));

% vFbp(PEP,FBP)
ratio = 1+M(fdp)/KFbpFBP;
flux(3) = vFbpmax.*(ratio-1).*(ratio).^3/(ratio.^4+Lfbp*(1+M(pep)./KFbpPEP).^(-4));

% vEX(PEP)
flux(4) = vEXmax.*M(pep)./(M(pep)+KEXPEP);

% vPEPout
flux(5) = kPEPout*M(pep);



