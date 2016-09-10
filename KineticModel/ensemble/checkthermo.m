function check = checkthermo(model,pvec,mc,Vind,fhandle)
% check thermodynamic consistency of given parameters
flux = zeros(model.nt_rxn,1);
check = zeros(model.nt_rxn,1);

flux(Vind) = fhandle(model,pvec,mc,Vind);
flxdelG = flux.*pvec.delGr;

nnzdelG = flxdelG(pvec.delGr~=0);
nnzcheck = zeros(length(nnzdelG),1);
nnzcheck(nnzdelG<0) = 1;
nnzcheck(nnzdelG>=0) = -1;
check(pvec.delGr~=0) = nnzcheck;

zdelG = flxdelG(pvec.delGr==0);
zcheck = zeros(length(zdelG),1);
zcheck(zdelG<1e-6) = 1;
zcheck(zdelG>=1e-6) = -1;
check(pvec.delGr==0) = zcheck;