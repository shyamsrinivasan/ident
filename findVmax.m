function pvec = findVmax(model,pvec,mc)
bmr = model.bmrxn;
Vuptake = model.Vuptake;
S = -model.S(:,bmr);
mu = 0.1;

%indices - fluxes
vglc = strcmpi(model.rxns,'exGLC');
vpts = strcmpi(model.rxns,'glcpts');
vpgi = strcmpi(model.rxns,'pgi');
vg6pd = strcmpi(model.rxns,'g6pdh2r');
vpgl = strcmpi(model.rxns,'pgl');
vgnd = strcmpi(model.rxns,'gnd');
vtkt1 = strcmpi(model.rxns,'tkt1');
vtkt2 = strcmpi(model.rxns,'tkt2');
vtala = strcmpi(model.rxns,'tala');
vrpi = strcmpi(model.rxns,'rpi');
vrpe = strcmpi(model.rxns,'rpe');
vpfk = strcmpi(model.rxns,'pfk');
vfba = strcmpi(model.rxns,'fba');
vtpi = strcmpi(model.rxns,'tpi');
vgapd = strcmpi(model.rxns,'gapd');

g6p = strcmpi(model.mets,'g6p[c]');
pgl = strcmpi(model.mets,'6pgl[c]');
pgc = strcmpi(model.mets,'6pgc[c]');
r5p = strcmpi(model.mets,'r5p[c]');
ru5p = strcmpi(model.mets,'ru5p-D[c]');
x5p = strcmpi(model.mets,'xu5p-D[c]');
s7p = strcmpi(model.mets,'s7p[c]');
e4p = strcmpi(model.mets,'e4p[c]');
f6p = strcmpi(model.mets,'f6p[c]');
fdp = strcmpi(model.mets,'fdp[c]');
dhap = strcmpi(model.mets,'dhap[c]');
g3p = strcmpi(model.mets,'g3p[c]');

%known fluxes
% knflux = {'G6PDH2r','me1','mdh'};
% knflxid = cellfun(@(x)strcmpi(model.rxns,x),knflux,'UniformOutput',false);
% knflxid = cell2mat(cellfun(@find,knflxid,'UniformOutput',false));
knflxid = logical(pvec.Vmax==0);
knid = false(model.nt_rxn,1);
knid(knflxid) = 1;

%initial flux calculation
vf = Vuptake;
vf(knid) = CKinetics(model,pvec,mc,find(knid));

vf(logical(Vuptake)) = Vuptake(logical(Vuptake));
knid(logical(Vuptake)) = 1;
vf(bmr) = mu;
knid(bmr) = 1;

[~,ck] = CKinetics(model,pvec,mc,find(vpts));
pvec.Vmax(vpts) = vf(vglc)/ck;
knid(vpts) = 1;
vf(vpts) = CKinetics(model,pvec,mc,find(vpts));

%upper glycolysis I   
[vf,pvec,knid] = linearfluxes(model,vf,mc,mu,knid,pvec,g6p,vpgi);

%pentose phosphate
nr = vf(vg6pd)-mu*mc(pgl);
[pvec,vf,knid] = Vmax_sub1(model,pvec,mc,find(vpgl),nr,vf,knid);

nr = vf(vpgl)-mu*mc(pgc);
[pvec,vf,knid] = Vmax_sub1(model,pvec,mc,find(vgnd),nr,vf,knid);

nr = (vf(vgnd)+mu*(-2*mc(s7p)-mc(r5p)-S(r5p)-mc(ru5p)+mc(e4p)+S(e4p)-mc(x5p)))/3;
[pvec,vf,knid] = Vmax_sub1(model,pvec,mc,find(vtala),nr,vf,knid);

nr = -vf(vtala)-mu*mc(s7p);
[pvec,vf,knid] = Vmax_sub1(model,pvec,mc,find(vtkt1),nr,vf,knid);

nr = -vf(vtala)+mu*mc(e4p)+mu*S(e4p);
[pvec,vf,knid] = Vmax_sub1(model,pvec,mc,find(vtkt2),nr,vf,knid);

nr = vf(vtkt1)-mu*mc(r5p)-mu*S(r5p);
[pvec,vf,knid] = Vmax_sub1(model,pvec,mc,find(vrpi),nr,vf,knid);

nr = vf(vrpi)+vf(vgnd)-mu*mc(ru5p);
[pvec,vf,knid] = Vmax_sub1(model,pvec,mc,find(vrpe),nr,vf,knid);

%upper glycolysis II
[vf,pvec,knid] = linearfluxes(model,vf,mc,mu,knid,pvec,f6p,vpfk);

[vf,pvec,knid] = linearfluxes(model,vf,mc,mu,knid,pvec,fdp,vfba);

[vf,pvec,knid] = linearfluxes(model,vf,mc,mu,knid,pvec,dhap,vtpi);

[vf,pvec,knid] = linearfluxes(model,vf,mc,mu,knid,pvec,g3p,vgapd);

%lower glycolysis 
[vf,pvec,knid] = linearfluxes(model,vf,mc,mu,knid,pvec,dpg,vpgk);

[vf,pvec,knid] = linearfluxes(model,vf,mc,mu,knid,pvec,pg3,vpgm);

[vf,pvec,knid] = linearfluxes(model,vf,mc,mu,knid,pvec,pg2,veno);

%how to calculate v(pyk) from assumed Vmax?
[vf,pvec,knid] = linearfluxes(model,vf,mc,mu,knid,pvec,pep,vppc);



function [vf,pvec,knid] = linearfluxes(model,vf,mc,mu,knid,pvec,metid,vid)
    vmet_f = model.S(metid,:);
    vknw = setdiff(find(vmet_f),find(vid));
    nr = model.S(metid,vknw)*vf(vknw)-mu*mc(metid);
    [pvec,vf,knid] = Vmax_sub1(model,pvec,mc,find(vid),nr,vf,knid);
%     [~,ck] = CKinetics(model,pvec,mc,find(rxnid));
%     pvec.Vmax(rxnid) = (model.S(metid,vknw)*vf(vknw)-mu*mc(metid))/ck;
%     knidx(rxnid) = 1;
%     vf(rxnid) = CKinetics(model,pvec,mc,find(rxnid));
return

function [pvec,vf,knid] = Vmax_sub1(model,pvec,mc,vid,nr,vf,knid)

[~,ck] = CKinetics(model,pvec,mc,vid);
if ck
    if nr/ck < 0
        pvec.Vmax(vid) = model.Vss(vid)/ck;
    else
        pvec.Vmax(vid) = nr/ck;
    end
else
    pvec.Vmax(vid) = 0;
end
vf(vid) = CKinetics(model,pvec,mc,vid);
knid(vid) = 1;