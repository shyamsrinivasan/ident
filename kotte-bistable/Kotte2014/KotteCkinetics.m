function dM = KotteCkinetics(kmrgd,pvec,model)
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
d = pvec(13);

allmc = [kmrgd;model.PM];

pep = strcmpi(model.mets,'pep[c]');
fdp = strcmpi(model.mets,'fdp[c]');
enz = strcmpi(model.mets,'enz[c]');
ac = strcmpi(model.mets,'ac[e]');

% convenience kinetics model
nr1 = allmc(enz)*kEcat*prod(allmc(ac)./KEacetate);
dr1 = 1 + prod(allmc(ac)./KEacetate);
act1 = 1;
flux(1) = nr1/dr1*act1; % J

nr2 = vemax;
dr2 = 1;
act2 = prod(1./(1+(allmc(fdp)./KeFBP).^ne));
flux(2) = nr2/dr2*act2; % E flux

nr3 = vFbpmax*prod(allmc(fdp)./KFbpFBP);
dr3 = 1 + prod(allmc(fdp)./KFbpFBP);
act3 = prod(1 + (allmc(pep)./KFbpPEP).^4);
flux(3) = nr3/dr3*act3; % FBP

nr4 = vEXmax*prod(allmc(pep)./KEXPEP);
dr4 = 1 + prod(allmc(pep)./KEXPEP);
act4 = 1;
flux(4) = nr4/dr4*act4; % X

% original model
% v(1) = kEcat.*allmc(enz).*acetate./(acetate+KEacetate);
% v(2) = vemax.*(1-1./(1+(KeFBP./allmc(fdp)).^ne));
% ratio = 1+allmc(fdp)/KFbpFBP;
% v(3) = vFbpmax.*(ratio-1).*(ratio).^3/(ratio.^4+Lfbp*(1+allmc(pep)./KFbpPEP).^(-4));
% v(4) = vEXmax.*allmc(pep)./(allmc(pep)+KEXPEP);

R(1) = v(1) - flux(1);
R(2) = v(2) - flux(2);
R(3) = v(3) - flux(3);
R(4) = v(4) - flux(4);

% dM = zeros(length(kmrgd),1);
% dM(1) = flux(1)-flux(4);    % pep
% dM(2) = flux(4)-flux(3);    % fdp
% dM(3) = flux(2)-d*allmc(3); % E

dM = R;
% dM(pep) = model.S(pep,:)*flux';
% dM(fdp) = model.S(fdp,:)*flux';
% dM(enz) = flux(tfr) - d*allmc(enz);

% if any(kmrgd<0)
%     status = -1;
% else
    status = 0;
% end

