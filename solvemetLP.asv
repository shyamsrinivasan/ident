function solvemetLP(model)
%
S = model.S;
Keq = model.Keq;
Vss = model.Vss;
rxns = model.rxns;
Vind = model.Vind;
Vex = model.Vex;
VFex = model.VFex;

vflux = zeros(size(model.S,2),1);
vflux([Vind Vex VFex']) = 1;

%remove reactions with zero fluxes
S(:,abs(Vss)<1e-6) = [];
Keq(abs(Vss)<1e-6) = [];
rxns(abs(Vss)<1e-6) = [];
Vss(abs(Vss)<1e-6) = [];
vflux(abs(Vss)<1e-6) = [];


%nonzero reactions for thermodynamic analysis
[nmet,nrxn] = size(S);

%Ax <=b 
A = S';
A = repmat(sign(Vss),1,nmet).*A;

%bounds
lb = zeros(nmet,1);
lb(lb==0) = log(1e-28);
ub = zeros(nmet,1);
ub(ub==0) = log(1e3);

%concentrations in M = mole/L, Bennett et al., 2009
lb(strcmpi(model.mets,'glc[e]')) = log(0.2);
lb(strcmpi(model.mets,'atp[c]')) = log(8.13e-3);
lb(strcmpi(model.mets,'adp[c]')) = log(4.37e-4);
lb(strcmpi(model.mets,'amp[c]')) = log(2.32e-4);
lb(strcmpi(model.mets,'pep[c]')) = log(1.46e-4);
lb(strcmpi(model.mets,'dhap[c]')) = log(3.44e-4);
lb(strcmpi(model.mets,'nad[c]')) = log(2.32e-3);
lb(strcmpi(model.mets,'nadh[c]')) = log(5.45e-5);
lb(strcmpi(model.mets,'nadp[c]')) = log(1.4e-7);
lb(strcmpi(model.mets,'nadph[c]')) = log(1.1e-4);
lb(strcmpi(model.mets,'accoa[c]')) = log(5.29e-4);
lb(strcmpi(model.mets,'coa[c]')) = log(8.8e-5);
lb(strcmpi(model.mets,'cit[c]')) = log(1.10e-3);
lb(strcmpi(model.mets,'mal[c]')) = log(1.66e-3);
lb(strcmpi(model.mets,'fum[c]')) = log(3e-6);
lb(strcmpi(model.mets,'succ[c]')) = log(3.41e-4);
lb(strcmpi(model.mets,'succoa[c]')) = log(1.42e-4);

ub(strcmpi(model.mets,'glc[e]')) = log(0.2);
ub(strcmpi(model.mets,'atp[c]')) = log(8.13e-3);
ub(strcmpi(model.mets,'adp[c]')) = log(4.37e-4);
ub(strcmpi(model.mets,'amp[c]')) = log(2.32e-4);
ub(strcmpi(model.mets,'pep[c]')) = log(1.46e-4);
ub(strcmpi(model.mets,'dhap[c]')) = log(3.44e-4);
ub(strcmpi(model.mets,'nad[c]')) = log(2.32e-3);
ub(strcmpi(model.mets,'nadh[c]')) = log(5.45e-5);
ub(strcmpi(model.mets,'nadp[c]')) = log(1.4e-7);
ub(strcmpi(model.mets,'nadph[c]')) = log(1.1e-4);
ub(strcmpi(model.mets,'accoa[c]')) = log(5.29e-4);
ub(strcmpi(model.mets,'coa[c]')) = log(8.8e-5);
ub(strcmpi(model.mets,'cit[c]')) = log(1.10e-3);
ub(strcmpi(model.mets,'mal[c]')) = log(1.66e-3);
ub(strcmpi(model.mets,'fum[c]')) = log(3e-6);
ub(strcmpi(model.mets,'succ[c]')) = log(3.41e-4);
ub(strcmpi(model.mets,'succoa[c]')) = log(1.42e-4);

b = sign(Vss).*log(Keq);

cprod = sparse(1,nmet);

[x,xobj,flag] = cplexlp(cprod(:),A,b,[],[],lb,ub);
