function flux = noisyflux_kotte(x,pstruct,flux)
if isfield(pstruct,'p')
    p = pstruct.p;
end
if isfield(pstruct,'model')
    model = pstruct.model;
end
% perturbation - deletion
if isfield(pstruct,'del')
    del = pstruct.del;
else
    del = [];
end
% adatpted from Kotte_givenFlux
% May 2017
if nargin < 3
    flux = zeros(5,size(x,2));
%     flux = cons(flux,x);
end
% if nargin < 3
%     model = struct([]);
% end

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
%     acetate = 2;        % a.u acetate
    kPEPout = 0.2;
else
    % parameters    
    KEacetate = p(1);    % or 0.02
    KFbpFBP = p(2);
    Lfbp = p(3);
    KFbpPEP = p(4);
    KEXPEP = p(5);
    vemax = p(6);        % for bifurcation analysis: 0.7:0.1:1.3
    KeFBP = p(7);        % or 0.45
    ne = p(8);             % or 2
    kPEPout = p(11);
    kEcat = p(12);   
    vFbpmax = p(13);    
    vEXmax = p(14);   
end

pep = strcmpi(model.mets,'pep[c]');
fdp = strcmpi(model.mets,'fdp[c]');
enz = strcmpi(model.mets,'enz[c]');
ac = strcmpi(model.mets,'ac[e]');

% metabolic fluxes
% J(E, acetate)
flux(1,:) = kEcat.*x(enz,:).*x(ac,:)./(x(ac,:)+KEacetate);

% enzyme production fluxes
% E(FBP) for J (%FBP ---| Cra and Cra ---> E)
flux(2,:) = vemax.*(1-1./(1+(KeFBP./x(fdp,:)).^ne));

% vFbp(PEP,FBP)
ratio = 1+x(fdp,:)./KFbpFBP;
flux(3,:) = vFbpmax.*(ratio-1).*(ratio).^3./...
            (ratio.^4+Lfbp.*(1+x(pep,:)./KFbpPEP).^(-4));

% vEX(PEP)
flux(4,:) = vEXmax.*x(pep,:)./(x(pep,:)+KEXPEP);

% vPEPout
flux(5,:) = kPEPout*x(pep,:);

% generate noise
% noise = random('Normal',0,sqrt(var(flux)),5,1);
noise = rand(5,1)*2;

flux = flux + repmat(noise,1,size(flux,2));
flux(del,:) = 0;

