function flux = Kotte_glycolysisflux(M,pvec,flux,model)
if nargin < 4
    model = struct([]);
end
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
    % parameters
    if ~isstruct(pvec)
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
if ~isstruct(pvec) || isempty(model)
    flux(1) = kEcat.*M(1).*acetate./(acetate+KEacetate);
end

% vEX(PEP)
if ~isstruct(pvec) || isempty(model)
    flux(4) = vEXmax.*M(2)./(M(2)+KEXPEP);
end

% % vFbp(PEP,FBP)

if ~isstruct(pvec) || isempty(model)
    ratio = 1+M(3)/KFbpFBP;
    flux(3) = vFbpmax.*(ratio-1).*(ratio).^3/(ratio.^4+Lfbp*(1+M(2)./KFbpPEP).^(-4));
else
%     tfr = strcmpi(model.rxns,'FBP'); 
%     tfm = strcmpi(model.mets,'fdp[c]');
%     tfrg = strcmpi(model.mets,'pep[c]');
%     dr_sb = prod((1+M(tfm)./pvec.K(tfm,tfr)).^1);    
%     
%     flux(tfr) = (pvec.kcat_fwd(tfr)*prod((M(tfm)./pvec.K(tfm,tfr)).^1).*...
%                 prod((M(tfrg)./pvec.KIact(tfrg,tfr)).^4./...
%                           ((1+M(tfrg)./pvec.KIact(tfrg,tfr)).^4)))/(dr_sb);                
end

% enzyme production fluxes
% E(FBP) for J (%FBP ---| Cra and Cra ---> E)
if ~isempty(model)
    ne = 2;
    tfr = strcmpi(model.rxns,'ENZC'); 
    tfrg = strcmpi(model.mets,'fdp[c]');
    flux(tfr) = pvec.Vmax(tfr).*(1-1./(1+pvec.KIact(tfrg,tfr)./M(tfrg).^ne));
else
    flux(2) = vemax.*(1-1./(1+(KeFBP./M(3)).^ne));
end