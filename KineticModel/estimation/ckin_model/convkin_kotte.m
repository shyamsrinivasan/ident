function flux = convkin_kotte(x,pstruct,flux)
if nargin < 3
    flux = zeros(5,size(x,2));
end
if isfield(pstruct,'p')
    p = pstruct.p;
end
if isfield(pstruct,'model')
    model = pstruct.model;
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
    K1pep = p(15);  % .02
    KEXfdp = p(16); % .3
    rhoA = p(17); % 0.5
end

pep = strcmpi(model.mets,'pep[c]');
fdp = strcmpi(model.mets,'fdp[c]');
enz = strcmpi(model.mets,'enz[c]');
ac = strcmpi(model.mets,'ac[e]');

% convinience kientics for acetate uptake and conversion to pep
flux(1,:) = kEcat.*x(enz,:).*(x(ac,:)./KEacetate)/...
            (1 + x(ac,:)./KEacetate + x(pep,:)./K1pep);

% hill kinetics for enzyme production        
flux(2,:) = vemax.*(1-1./(1+(KeFBP./x(fdp,:)).^ne));

% convinience kinetics w/ allostery for vFbp 
acratio = x(pep,:)./KFbpPEP;
acflx = rhoA + (1-rhoA).*(acratio./(1+acratio)).^4;
% acflx = (rhoA + (1-rhoA).*acratio./(1+acratio)).^4;
flux(3,:) = vFbpmax.*acflx.*(x(fdp,:)./KFbpFBP)./(1+x(fdp,:)./KFbpFBP);

% convinience kinetics for vEX(pep)
flux(4,:) = vEXmax.*(x(pep,:)./KEXPEP)/(1+x(pep,:)./KEXPEP+x(fdp,:)./KEXfdp);

% linear kinetics for vPEPout
flux(5,:) = kPEPout*x(pep,:);











