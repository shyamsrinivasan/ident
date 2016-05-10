function J = KottegivenJacobian(M,pvec,model)

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

allmc = [M;model.PM];
    
flux = Kotte_givenFlux(allmc,pvec,model);
pep = strcmpi(model.mets,'pep[c]');
fdp = strcmpi(model.mets,'fdp[c]');
enz = strcmpi(model.mets,'enz[c]');
ac = strcmpi(model.mets,'ac[e]');

J = sparse(3,3);

df4dx1 = flux(4)*(1/allmc(pep)-1/(allmc(pep)+KEXPEP));
df1dx3 = flux(1)*allmc(enz);
ratio = 1+allmc(fdp)/KFbpFBP;
Drf3 = ratio.^4+Lfbp*(1+allmc(pep)./KFbpPEP).^(-4);
df3dx1 = flux(3)*(Lfbp/Drf3)*(1+allmc(pep)/KFbpPEP)^(-5)*(4/KFbpPEP);
df2dx2 = -2*flux(2)/allmc(fdp)*(1-flux(2)/vemax);
df3dx2 = flux(3)*(1/allmc(fdp)+3/(allmc(fdp)+KFbpFBP)-4*flux(3)/vFbpmax*(1/allmc(fdp))/(1+allmc(fdp)/KFbpFBP));

J(1,1) = -df4dx1-0.2;
J(1,3) = df1dx3;
J(2,1) = df4dx1-df3dx1;
J(2,2) = -df3dx2;
J(3,2) = df2dx2;



