function x = initialsample(model)
%
S = model.S;
Keq = model.Keq;
Vss = model.Vss;

%remove reactions with zero fluxes
S(:,abs(Vss)<1e-6) = [];
Keq(abs(Vss)<1e-6) = [];
Vss(abs(Vss)<1e-6) = [];

%nonzero reactions for thermodynamic analysis
[nmet,nrxn] = size(S);

%Ax <=b 
A = S';
vsign = sign(Vss);
A = repmat(sign(Vss),1,nmet).*A;

x = zeros(nmet,1);
b = vsign.*log(Keq);
lb = zeros(nmet,1);
lb(lb==0) = log(1e-6);
ub = zeros(nmet,1);
ub(ub==0) = log(1000);
cprod = sparse(1,nmet);

[x,xobj,flag] = cplexlp(cprod(:),A,b,[],[],lb,ub);

if flag > 0
    x = x;
end
