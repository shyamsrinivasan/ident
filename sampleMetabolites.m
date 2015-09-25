function mc = sampleMetabolites(model,mc)
if nargin < 2
    mc = zeros(model.nt_metab,1);
end

delG = model.delSGr;
R = 0.008314; %kJ/mol.K
T = 298; %K
RT = R*T;

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
vtala = strcmpi(model.rxns,'tala');
vtkt2 = strcmpi(model.rxns,'tkt2');
vpfk = strcmpi(model.rxns,'pfk');
% vfbp = strcmpi(model.rxns,'fbp');
vfba = strcmpi(model.rxns,'fba');
vtpi = strcmpi(model.rxns,'tpi');
vgapd = strcmpi(model.rxns,'gapd');
vpgk = strcmpi(model.rxns,'pgk');
vpgm = strcmpi(model.rxns,'pgm');
% veno = strcmpi(model.rxns,'eno');
vacont = strcmpi(model.rxns,'aconta');
vcs = strcmpi(model.rxns,'cs');
vicl = strcmpi(model.rxns,'icl');
vicd = strcmpi(model.rxns,'icdhyr');
% vakgd = strcmpi(model.rxns,'akgdh');
% vsuca = strcmpi(model.rxns,'sucoas');
% vsucd = strcmpi(model.rxns,'sucdi');
% vfum = strcmpi(model.rxns,'fum');
vmals = strcmpi(model.rxns,'mals');
% vme1 = strcmpi(model.rxns,'me1');
vmdh = strcmpi(model.rxns,'mdh');
% vpdh = strcmpi(model.rxns,'pdh');
vpfl = strcmpi(model.rxns,'pfl');
% vppc = strcmpi(model.rxns,'ppc');
% vppck = strcmpi(model.rxns,'ppck');
vpyk = strcmpi(model.rxns,'pyk');
vldh = strcmpi(model.rxns,'ldh_d');
vpta = strcmpi(model.rxns,'ptar');
vack = strcmpi(model.rxns,'ackr');
vacald = strcmpi(model.rxns,'acald');
valcd2x = strcmpi(model.rxns,'alcd2x');

%glycolysis
mcub(pyr)=mc(pep)*mc(adp)/mc(atp)*exp(-delG(vpyk)/RT);
mc(pyr) = random(makedist('Uniform','lower',mcub(pyr)/10,'upper',mcub(pyr)),1,1);

mcub(g6p)=mc(glc)*mc(pep)/mc(pyr)*exp(-delG(vpts)/RT);
mc(g6p) = random(makedist('Uniform','lower',mcub(g6p)/100,'upper',mcub(g6p)/10),1,1);

mcub(f6p) = mc(g6p)*exp(-delG(vpgi)/RT);
mc(f6p) = random(makedist('Uniform','lower',mcub(f6p)/10,'upper',mcub(f6p)),1,1);

mcub(fdp) = mc(f6p)*exp(-delG(vpfk)/RT);
mc(fdp) = random(makedist('Uniform','lower',mcub(fdp)/10,'upper',mcub(fdp)),1,1);

mcub(g3p) = min(mc(fdp)/mc(dhap)*exp(-delG(vfba)/RT),...
                mc(dhap)*exp(-delG(vtpi)/RT));
mc(g3p) = random(makedist('Uniform','lower',mcub(g3p)/10,'upper',mcub(g3p)),1,1);

mcub(dpg) = mc(g3p)*mc(nad)/mc(nadh)*exp(-delG(vgapd)/RT);
mc(dpg) = random(makedist('Uniform','lower',mcub(dpg)/10,'upper',mcub(dpg)),1,1);

mcub(pg3) = mc(dpg)*mc(adp)/mc(atp)*exp(-delG(vpgk)/RT);
mc(pg3) = random(makedist('Uniform','lower',mcub(pg3)/10,'upper',mcub(pg3)),1,1);

mcub(pg2) = mc(pg3)*exp(-delG(vpgm)/RT);
mc(pg2) = random(makedist('Uniform','lower',mcub(pg2)/10,'upper',mcub(pg2)),1,1);

%pentose phosphate
mc(pgl) = mc(g6p)*mc(nadp)/mc(nadph)*exp(-delG(vg6pd)/RT);
%mc(pgl)  = random(makedist('Uniform','lower',mcub(pgl)/10,'upper',mcub(pgl),1,1));

mc(pgc) = mc(pgl)*exp(-delG(vpgl)/RT);
% mc(pgc) = random(makedist('Uniform','lower',mcub(pgc)/10,'upper',mcub(pgc),1,1));

mc(ru5p) = mc(pgc)*mc(nadp)/mc(nadph)*exp(-delG(vgnd)/RT);
% mc(ru5p) = random(makedist('Uniform','lower',mcub(ru5p)/10,'upper',mcub(ru5p),1,1));

mclb(x5p) = mc(ru5p)*exp(-delG(vrpe)/RT);
mc(x5p) = random(makedist('Uniform','lower',mclb(x5p),'upper',mclb(x5p)*10),1,1);

mclb(r5p) = mc(ru5p)*exp(-delG(vrpi)/RT);
mc(r5p) = random(makedist('Uniform','lower',mclb(r5p),'upper',mclb(r5p)*10),1,1);

mclb(s7p) = mc(r5p)*mc(x5p)/mc(g3p)*exp(-delG(vtkt1)/RT);
mc(s7p) = random(makedist('Uniform','lower',mclb(s7p),'upper',mclb(s7p)*10),1,1);

mcub(e4p) = mc(f6p)*mc(g3p)/mc(x5p)*exp(delG(vtkt2)/RT);
mclb(e4p) = mc(g3p)*mc(s7p)/mc(f6p)*exp(-delG(vtala)/RT);
mc(e4p) = random(makedist('Uniform','lower',mclb(e4p),'upper',mcub(e4p)),1,1);

%tca cycle
mclb(oaa) = mc(cit)*mc(coa)/mc(accoa)*exp(delG(vcs)/RT);
mcub(oaa) = mc(mal)*mc(nad)/mc(nadh)*exp(delG(vmdh)/RT);
mc(oaa) = random(makedist('Uniform','lower',mclb(oaa),'upper',mcub(oaa)),1,1);

mcub(icit) = mc(cit)*exp(-delG(vacont)/RT);
mc(icit) = random(makedist('Uniform','lower',mcub(icit)/10,'upper',mcub(icit)),1,1);

mcub(akg) = mc(icit)*mc(nadp)/mc(nadph)*exp(-delG(vicd)/RT);
mc(akg) = random(makedist('Uniform','lower',mcub(akg)/10,'upper',mcub(akg)),1,1);

%glyoxylate
mcub(glx) = mc(icit)/mc(suc)*exp(-delG(vicl)/RT);
mclb(glx) = mc(coa)*mc(mal)/mc(accoa)*exp(delG(vmals)/RT);
mc(glx) = random(makedist('Uniform','lower',mclb(glx),'upper',mcub(glx)),1,1);

%pyruvate metabolism
mcub(form) = mc(coa)*mc(pyr)/mc(accoa)*exp(-delG(vpfl)/RT);
mc(form) = random(makedist('Uniform','lower',mcub(form)/10,'upper',mcub(form)),1,1);

mcub(lac) = mc(pyr)*mc(nadh)/mc(nad)*exp(-delG(vldh)/RT);
mc(lac) = random(makedist('Uniform','lower',mcub(lac)/10,'upper',mcub(lac)),1,1);

mcub(actp) = mc(accoa)/mc(coa)*exp(-delG(vpta)/RT);
mc(actp) = random(makedist('Uniform','lower',mcub(actp)/10,'upper',mcub(actp)),1,1);

mcub(ac) = mc(actp)*mc(adp)/mc(atp)*exp(delG(vack)/RT);
mc(ac) = random(makedist('Uniform','lower',mcub(ac)/10,'upper',mcub(ac)),1,1);

mcub(acald) = mc(accoa)*mc(nadh)/(mc(coa)*mc(nad))*exp(delG(vacald)/RT);
mc(acald) = random(makedist('Uniform','lower',mcub(acald)/10,'upper',mcub(acald)),1,1);

mcub(etoh) = mc(acald)*mc(nadh)/mc(nad)*exp(-delG(valcd2x)/RT);
mc(etoh) = random(makedist('Uniform','lower',mcub(etoh)/10,'upper',mcub(etoh)),1,1);

%convert M to mM
mc = mc*1e3;













