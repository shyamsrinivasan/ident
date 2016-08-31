function flux = Kotte_CFlux(M,pvec,model,flux)
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
end

pep = strcmpi(model.mets,'pep[c]');
fdp = strcmpi(model.mets,'fdp[c]');
enz = strcmpi(model.mets,'enz[c]');
ac = strcmpi(model.mets,'ac[e]');

% convenience kinetics model
nr1 = M(enz)*kEcat*prod(M(ac)./KEacetate);
dr1 = 1 + prod(M(ac)./KEacetate);
act1 = 1;
flux(1) = nr1/dr1*act1; % J

nr2 = vemax;
dr2 = 1;
act2 = prod(1./(1+(M(fdp)./KeFBP).^ne));
flux(2) = nr2/dr2*act2; % E flux

nr3 = vFbpmax*prod(M(fdp)./KFbpFBP);
dr3 = 1 + prod(M(fdp)./KFbpFBP);
act3 = prod(1 + (M(pep)./KFbpPEP).^4);
flux(3) = nr3/dr3*act3; % FBP

nr4 = vEXmax*prod(M(pep)./KEXPEP);
dr4 = 1 + prod(M(pep)./KEXPEP);
act4 = 1;
flux(4) = nr4/dr4*act4; % vEX(PEP)