function [mc,delGr] = metsample(model,mc)
if nargin < 2
    mc = zeros(model.nt_metab,1);
end

delSG = model.delSGr;
delGr = delGsample(model);

R = 0.008314; %kJ/mol.K
T = 298; %K
RT = R*T;

Q = exp((delGr-delSG)./RT);

%metabolites
glc = strcmpi(model.mets,'glc[e]');
g6p = strcmpi(model.mets,'g6p[c]');
f6p = strcmpi(model.mets,'f6p[c]');
pgl = strcmpi(model.mets,'6pgl[c]');
pgc = strcmpi(model.mets,'6pgc[c]');
r5p = strcmpi(model.mets,'r5p[c]');
ru5p = strcmpi(model.mets,'ru5p-D[c]');
fdp = strcmpi(model.mets,'fdp[c]');
dhap = strcmpi(model.mets,'dhap[c]');
g3p = strcmpi(model.mets,'g3p[c]');
dpg = strcmpi(model.mets,'13dpg[c]');
x5p = strcmpi(model.mets,'xu5p-D[c]');
e4p = strcmpi(model.mets,'e4p[c]');
s7p = strcmpi(model.mets,'s7p[c]');
pg3 = strcmpi(model.mets,'3pg[c]');
pg2 = strcmpi(model.mets,'2pg[c]');
pep = strcmpi(model.mets,'pep[c]');
pyr = strcmpi(model.mets,'pyr[c]');
accoa = strcmpi(model.mets,'accoa[c]');
coa = strcmpi(model.mets,'coa[c]');
cit = strcmpi(model.mets,'cit[c]');
icit = strcmpi(model.mets,'icit[c]');
glx = strcmpi(model.mets,'glx[c]');
akg = strcmpi(model.mets,'akg[c]');
sucoa = strcmpi(model.mets,'succoa[c]');
suc = strcmpi(model.mets,'succ[c]');
fum = strcmpi(model.mets,'fum[c]');
mal = strcmpi(model.mets,'mal[c]');
oaa = strcmpi(model.mets,'oaa[c]');
form = strcmpi(model.mets,'for[c]');
lac = strcmpi(model.mets,'lac[c]');
actp = strcmpi(model.mets,'actp[c]');
ac = strcmpi(model.mets,'ac[c]');
acald = strcmpi(model.mets,'acald[c]');
etoh = strcmpi(model.mets,'etoh[c]');

nad = strcmpi(model.mets,'nad[c]');
nadh = strcmpi(model.mets,'nadh[c]');
nadp = strcmpi(model.mets,'nadp[c]');
nadph = strcmpi(model.mets,'nadph[c]');
atp = strcmpi(model.mets,'atp[c]');
adp = strcmpi(model.mets,'adp[c]');
amp = strcmpi(model.mets,'amp[c]');
hc = strcmpi(model.mets,'h[c]');
he = strcmpi(model.mets,'h[e]');
q8 = strcmpi(model.mets,'q8[c]');
q8h2 = strcmpi(model.mets,'q8h2[c]');
pi = strcmpi(model.mets,'pi[c]');

%concentrations in M = mole/L, Bennett et al., 2009
mc(glc) = 0.2;
mc(atp) = 8.13e-3;
mc(adp) = 4.37e-4;
mc(amp) = 2.32e-4;
mc(pep) = 1.46e-4;
mc(dhap) = 3.44e-4;
mc(nad) = 2.32e-3;
mc(nadh) = 5.45e-5;
mc(nadp) = 1.4e-7;
mc(nadph) = 1.1e-4;
mc(accoa) = 5.29e-4;
mc(coa) = 8.8e-5;
mc(cit) = 1.10e-3;
mc(mal) = 1.66e-3;
mc(fum) = 3e-6;
mc(suc) = 3.41e-4;
mc(sucoa) = 1.42e-4;

%neutral pH
mc(hc) = 1e-7;
mc(he) = 1e-7;

%water - 1 L
mc(strcmpi(model.mets,'h2o[c]')) = 55.0;

%free phosphate
mc(pi) = mc(atp)-mc(adp)-mc(amp)/10;

%dissolved co2
mc(strcmpi(model.mets,'co2[c]')) = 1e-6;

%dissolved o2
mc(strcmpi(model.mets,'o2[c]')) = 1e-5;

%indices - fluxes
% vglc = strcmpi(model.rxns,'exGLC');
vpts = strcmpi(model.rxns,'glcpts');
vg6pd = strcmpi(model.rxns,'G6PDH2r');
vpgi = strcmpi(model.rxns,'pgi');
vpgl = strcmpi(model.rxns,'pgl');
vgnd = strcmpi(model.rxns,'gnd');
vrpe = strcmpi(model.rxns,'rpe');
vrpi = strcmpi(model.rxns,'rpi');
vtkt1 = strcmpi(model.rxns,'tkt1');
vtkt2 = strcmpi(model.rxns,'tkt2');
vpfk = strcmpi(model.rxns,'pfk');
vtpi = strcmpi(model.rxns,'tpi');
vgapd = strcmpi(model.rxns,'gapd');
vpgk = strcmpi(model.rxns,'pgk');
vpgm = strcmpi(model.rxns,'pgm');
vacont = strcmpi(model.rxns,'aconta');
vcs = strcmpi(model.rxns,'cs');
vicl = strcmpi(model.rxns,'icl');
vicd = strcmpi(model.rxns,'icdhyr');
vmals = strcmpi(model.rxns,'mals');
vpfl = strcmpi(model.rxns,'pfl');
vpyk = strcmpi(model.rxns,'pyk');
vldh = strcmpi(model.rxns,'ldh_d');
vpta = strcmpi(model.rxns,'ptar');
vack = strcmpi(model.rxns,'ackr');
vacald = strcmpi(model.rxns,'acald');
valcd2x = strcmpi(model.rxns,'alcd2x');
vsucd = strcmpi(model.rxns,'sucdi');

%glycolysis
mc(pyr)=mc(pep)*mc(adp)/(mc(atp)*Q(vpyk));
mc(g6p) = mc(glc)*mc(pep)/mc(pyr)*Q(vpts);
mc(f6p) = mc(g6p)*Q(vpgi);
mc(fdp) = mc(f6p)*mc(atp)/mc(adp)*Q(vpfk);
mc(g3p) = mc(dhap)*Q(vtpi);
mc(dpg) = mc(g3p)*mc(nad)*mc(pi)/mc(nadh)*Q(vgapd);
mc(pg3) = mc(dpg)*mc(adp)/(mc(atp)*Q(vpgk));
mc(pg2) =  mc(pg3)/Q(vpgm);

%pentose phosphate
mc(pgl) = mc(g6p)*mc(nadp)/mc(nadph)*Q(vg6pd);
mc(pgc) = mc(pgl)*Q(vpgl);
mc(ru5p) = mc(pgc)*mc(nadp)/mc(nadph)*Q(vgnd);
mc(x5p) = mc(ru5p)*Q(vrpe);
mc(r5p) = mc(ru5p)*Q(vrpi);
mc(s7p) = mc(r5p)*mc(x5p)/(mc(g3p)*Q(vtkt1));
mc(e4p) = mc(f6p)*mc(g3p)/mc(x5p)*Q(vtkt2);

%tca cycle
%ubiquinone/ubiquinol
if delGr(vsucd) ~= 0
    q8h2_q8 = mc(suc)/mc(fum)*Q(vsucd);
elseif delGr(vfrd) ~= 0
    q8h2_q8 = mc(suc)/mc(fum)*Q(vfrd7);
else
    q8h2_q8 = mc(suc)/mc(fum)*Q(vsucd);
end
mc(q8h2) = 1e-6;
mc(q8) = mc(q8h2)/q8h2_q8;

mc(oaa) = mc(cit)*mc(coa)/mc(accoa)*Q(vcs);
mc(icit) = mc(cit)*Q(vacont);
mc(akg) = mc(icit)*mc(nadp)/mc(nadph)*Q(vicd);

%glyoxylate
mc(glx) = mc(icit)/mc(suc)*Q(vicl);
mc(glx) = mc(coa)*mc(mal)/mc(accoa)*Q(vmals);

%pyruvate metabolism
mc(form) = mc(coa)*mc(pyr)/(mc(accoa)*Q(vpfl));
mc(lac) = mc(pyr)*mc(nadh)/(mc(nad)*Q(vldh));
mc(actp) = mc(accoa)*mc(pi)/mc(coa)*Q(vpta);
mc(ac) = mc(actp)*mc(adp)/(mc(atp)*Q(vack));
mc(acald) = mc(accoa)*mc(nadh)/(mc(coa)*mc(nad)*Q(vacald));
mc(etoh) = mc(acald)*mc(nadh)/(mc(nad)*Q(valcd2x));

%convert M to mM
mc = mc*1e3;








