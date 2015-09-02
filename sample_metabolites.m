function mc = sample_metabolites(model)
delG = model.delSGr;
R = 0.008314;%kJ/mol.K
T = 298;%K
RT = R*T;
%indices - fluxes
vglc = strcmpi(model.rxns,'exGLC');
vpts = strcmpi(model.rxns,'glcpts');
vg6pd = strcmpi(model.rxns,'G6PDH2r');
vpgi = strcmpi(model.rxns,'pgi');
vpgl = strcmpi(model.rxns,'pgl');
vgnd = strcmpi(model.rxns,'gnd');
vrpe = strcmpi(model.rxns,'rpe');
vfba = strcmpi(model.rxns,'fba');
vtpi = strcmpi(model.rxns,'tpi');
vgapd = strcmpi(model.rxns,'gapd');
vpgk = strcmpi(model.rxns,'pgk');
vpgm = strcmpi(model.rxns,'pgm');
veno = strcmpi(model.rxns,'eno');
vacont = strcmpi(model.rxns,'aconta');
vcs = strcmpi(model.rxns,'cs');
vicl = strcmpi(model.rxns,'icl');
vicd = strcmpi(model.rxns,'icdhyr');
vakgd = strcmpi(model.rxns,'akgdh');
vsuca = strcmpi(model.rxns,'sucoas');
vsucd = strcmpi(model.rxns,'sucdi');
vfum = strcmpi(model.rxns,'fum');
vmals = strcmpi(model.rxns,'mals');
vme1 = strcmpi(model.rxns,'me1');
vmdh = strcmpi(model.rxns,'mdh');
vpdh = strcmpi(model.rxns,'pdh');
vppc = strcmpi(model.rxns,'ppc');
vpyk = strcmpi(model.rxns,'pyk');
vatpm = strcmpi(model.rxns,'atpm');

%metabolites
glc = strcmpi(model.mets,'glc[e]');
g6p = strcmpi(model.mets,'g6p[c]');
pgl = strcmpi(model.mets,'6pgl[c]');
pgc = strcmpi(model.mets,'6pgc[c]');
r5p = strcmpi(model.mets,'r5p[c]');
ru5p = strcmpi(model.mets,'ru5p-D[c]');
fdp = strcmpi(model.mets,'fdp[c]');
dhap = strcmpi(model.mets,'dhap[c]');
g3p = strcmpi(model.mets,'g3p[c]');
dpg = strcmpi(model.mets,'13dpg[c]');
x5p = strcmpi(model.mets,'xu5p-D[c]');
s7p = strcmpi(model.mets,'s7p[c]');
pg3 = strcmpi(model.mets,'3pg[c]');
pg2 = strcmpi(model.mets,'2pg[c]');
pep = strcmpi(model.mets,'pep[c]');
pyr = strcmpi(model.mets,'pyr[c]');
accoa = strcmpi(model.mets,'accoa[c]');
cit = strcmpi(model.mets,'cit[c]');
icit = strcmpi(model.mets,'icit[c]');
glx = strcmpi(model.mets,'glx[c]');
akg = strcmpi(model.mets,'akg[c]');
sucoa = strcmpi(model.mets,'succoa[c]');
suc = strcmpi(model.mets,'succ[c]');
fum = strcmpi(model.mets,'fum[c]');
mal = strcmpi(model.mets,'mal[c]');
oaa = strcmpi(model.mets,'oaa[c]');
pic = strcmpi(model.mets,'pi[c]');
h = strcmpi(model.mets,'h[c]');

mc = zeros(model.nt_metab,1);
mclb = zeros(model.nt_metab,1);
mcub = zeros(model.nt_metab,1);
%assuymed concentrations in moles/L or M
mc(glc) = 0.2;
mc(pic) = 1e-3;
mc(h) = 1e-7;
nadp_r = 1e-3;%nadp/nadph
nad_r = 1e3;%nad/nadh

%g6p
mcub(g6p) = mc(glc)*mc(pi)*exp((-delG(vpts)+delG(vatpm)+delG(vpyk))/RT);
mc(g6p)
%f6p
mcub(f6p) = mc(g6p)*exp(-delG(vpgi)/RT);
mc(f6p)
%fdp
mclb(fdp) = mc(f6p)*mc(pic)*exp(delG(vfbp)/RT);
mcub(fdp) = mc(f6p)*mc(pic)*exp(delG(vfbp)/RT);
mc(fdp)
%atp/adp
atp_adplb = max(mc(h)*mc(pic)*exp(delG(vatpm)/RT),...
                mc(fdp)*mc(h)/mc(f6p)*exp(delG(vpfk)/RT));
%adp/atp
adp_atpub = min(1/(mc(h)*mc(pic))*exp(-delG(vatpm)/RT),...
                mc(f6p)/(mc(fdp)*mc(h))*exp(-delG(vpfk)/RT));  
adp_atp = 1/(mc(h)*mc(pic))*exp((-delG(vppc)-delG(vppck))/RT);            
%g3p
mcub(g3p) = sqrt(mc(fdp)*exp((-delG(vfba)-delG(-vtpi))/RT));
mc(g3p)
%dhap
mcub(dhap) = mc(g3p)*mc(fdp)*exp(-delG(vfba)/RT);
mclb(dhap) = mc(g3p)*exp(delG(vtpi)/RT);
mc(dhap)
%6pgl
mcub(pgl) = (mc(g6p)*nadp_r/mc(h))*exp(-delG(vg6pd)/RT);
mc(pgl)
%6pgc 
mcub(pgc) = (mc(pgl)/mc(h))*exp(-delG(vpgl)/RT);
mc(pgc)
%ru5p
mcub(ru5p) = mc(pgc)*nadp_r*exp(-delG(vgnd)/RT);
mc(ru5p)
%r5p
mcub(r5p) = mc(ru5p)*exp(delG(vrpi)/RT);
mc(r5p)
%xu5p
mcub(x5p) = mc(ru5p)*exp(-delG(vrpe)/RT);
mc(x5p)
%e4p
mclb(e4p) = (mc(f6p)*mc(g3p)/mc(x5p))*exp(delG(vtkt2)/RT);
mc(e4p)
%s7p
mcub(s7p) = (mc(x5p)*mc(r5p)/mc(g3p))*exp(-delG(vtkt1)/RT);
mclb(s7p) = (mc(e4p)*mc(f6p)/mc(g3p))*exp(delG(vtala)/RT);
mc(s7p)
%13dpg
mcub(dpg) = (mc(g3p)*mc(pic)*nad_r/mc(h))*exp(-delG(vgapd)/RT);
mc(dpg)
%3pg
mcub(pg3) = mc(dpg)*
mc(pg3)
%2pg
mcub(pg2) = mc(pg3)*exp(-delG(vpgm)/RT);
mc(pg2)
%pep
mcub(pep) = mc(pg2)*exp(-delG(eno)/RT);
mc(pep)
%pyr
mcub(pyr) = mc(pep)*mc(h)*exp(-delG(vpyk)/RT)*adp_atp;
mc(pyr)
%accoa



