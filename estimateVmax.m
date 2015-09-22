function pvec = estimateVmax(model,pvec,met)
bmrx = model.bmrxn;
Vuptake = model.Vuptake;
S = -model.S(:,bmrx);
%known fluxes
knflux = {'G6PDH2r','pdh','me1','mdh'};
knflx_id = cellfun(@(x)strcmpi(model.rxns,x),knflux,'UniformOutput',false);
knflx_id = cell2mat(cellfun(@find,knflx_id,'UniformOutput',false));
knidx = false(model.nt_rxn,1);
knidx(knflx_id) = 1;
%initial flux calculation
flux = Vuptake;
flux(knidx) = ConvinienceKinetics(model,pvec,met,find(knidx));
knidx(logical(Vuptake)) = 1;
n_known = length(find(knidx));
%indices - fluxes
vglc = strcmpi(model.rxns,'exGLC');
vpts = strcmpi(model.rxns,'glcpts');
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

%metabolites
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

%unknown fluxes
unkflux = {'glcpts','pgi','pgl','gnd','fba','rpe','rpi','tkt1',...
           'tala','tkt2','tpi','gapd','pgk','pgm','eno','ACONTa','cs',...
           'fum','sucdi','mals','icl','icdh','akgdh','SUCOAS',...
           'mdh','ppc','pyk'};  
%initial mu guess
mu_g = 0.5;
%initial mu calculated
mu_c = 0.1;%biomass_flux(model,met,[],Vuptake);
flux(bmrx) = mu_c;
Vmax = zeros(model.nt_rxn,1);

% while abs(mu_c-mu_g)>=1e-4
for i = 1:1
    %glc[e]    
    [~,vflux] = ConvinienceKinetics(model,pvec,met,find(vpts));
    pvec.Vmax(vpts) = flux(vglc)/vflux;
    knidx(vpts) = 1;
    flux(vpts) = ConvinienceKinetics(model,pvec,met,find(vpts));
    %g6p[c] - pgi       
    [flux,pvec,knidx] = linearfluxes(model,flux,met,mu_c,knidx,pvec,g6p,vpgi);
    %6pgl - pgl    
    [flux,pvec,knidx] = linearfluxes(model,flux,met,mu_c,knidx,pvec,pgl,vpgl);
    %6pgc - gnd    
    [flux,pvec,knidx] = linearfluxes(model,flux,met,mu_c,knidx,pvec,pgc,vgnd);
    %pfk,rpi,tkt1,tala,tkt2
    [flux,pvec,knidx] = pathway_pp(model,S,pvec,flux,met,mu_c,bmrx,knidx);
    %ru5p - vrpe
    [flux,pvec,knidx] = linearfluxes(model,flux,met,mu_c,knidx,pvec,ru5p,vrpe);    
    %fdp - vfba
    [flux,pvec,knidx] = linearfluxes(model,flux,met,mu_c,knidx,pvec,fdp,vfba);    
    %dhap - vtpi
    [flux,pvec,knidx] = linearfluxes(model,flux,met,mu_c,knidx,pvec,dhap,vtpi); 
    %g3p - vgapd
    [~,vflux] = ConvinienceKinetics(model,pvec,met,find(vgapd));
    pvec.Vmax(vgapd) = (flux(vfba)+flux(vgnd)+flux(vtpi)-mu_c*(S(g3p,1)+S(r5p,1)+...
                       met(g3p)+met(r5p)+met(ru5p)-met(s7p)+met(x5p)))/...
                       vflux;
    flux(vgapd) = ConvinienceKinetics(model,pvec,met,find(vgapd));
    knidx(vgapd) = 1;
    %13dpg - vpgk
    [flux,pvec,knidx] = linearfluxes(model,flux,met,mu_c,knidx,pvec,dpg,vpgk); 
    %3pg - vpgm
    [flux,pvec,knidx] = linearfluxes(model,flux,met,mu_c,knidx,pvec,pg3,vpgm); 
    %2pg - veno
    [flux,pvec,knidx] = linearfluxes(model,flux,met,mu_c,knidx,pvec,pg2,veno);     
    %vacont
    [~,vflux] = ConvinienceKinetics(model,pvec,met,find(vacont));
    pvec.Vmax(vacont) = (flux(veno)-mu_c*(met(pep)+met(pyr)-3*met(accoa)-...
             2*met(cit)+4*met(icit)+4*met(glx)+met(akg)+met(sucoa)+...
             met(fum)+met(mal)+met(oaa)))/vflux;
    flux(vacont) = ConvinienceKinetics(model,pvec,met,find(vacont));
    knidx(vacont) = 1;
    %vcs
    [~,vflux] = ConvinienceKinetics(model,pvec,met,find(vcs));
    pvec.Vmax(vcs) = (flux(vacont)+mu_c*met(cit))/vflux;
    flux(vcs) = ConvinienceKinetics(model,pvec,met,find(vcs));
    knidx(vcs) = 1;
    %vfum    
    [~,vflux] = ConvinienceKinetics(model,pvec,met,find(vfum));
    pvec.Vmax(vfum) = (flux(vacont)-mu_c*(met(icit)+met(akg)+met(sucoa)+...
                       met(suc)+met(fum)))/vflux;
    flux(vfum) = ConvinienceKinetics(model,pvec,met,find(vfum));
    knidx(vfum) = 1;
    %vsucdi
    [~,vflux] = ConvinienceKinetics(model,pvec,met,find(vsucd));
    pvec.Vmax(vsucd) = (flux(vfum)+mu_c*met(fum))/vflux;
    flux(vsucd) = ConvinienceKinetics(model,pvec,met,find(vsucd));
    knidx(vsucd) = 1;
    %vmals
    [~,vflux] = ConvinienceKinetics(model,pvec,met,find(vmals));
    pvec.Vmax(vmals) = (flux(vpdh)-flux(vcs)-mu_c*met(accoa))/vflux;
    flux(vmals) = ConvinienceKinetics(model,pvec,met,find(vmals));
    knidx(vmals) = 1;
    %vicl
    [~,vflux] = ConvinienceKinetics(model,pvec,met,find(vicl));
    pvec.Vmax(vicl) = (flux(vpdh)-flux(vcs)+mu_c*met(glx)-mu_c*met(accoa))/vflux;
    flux(vicl) = ConvinienceKinetics(model,pvec,met,find(vicl));
    knidx(vicl) = 1;
    %vicdh
    [~,vflux] = ConvinienceKinetics(model,pvec,met,find(vicd));
    pvec.Vmax(vicd) = (2*flux(vacont)-flux(vpdh)-...
                mu_c*(-met(glx)+mu_c*met(accoa)+met(cit)-met(icit)))/vflux;
    flux(vicd) = ConvinienceKinetics(model,pvec,met,find(vicd));
    knidx(vicd) = 1;
    %vakgd
    [~,vflux] = ConvinienceKinetics(model,pvec,met,find(vakgd));
    pvec.Vmax(vakgd) = (flux(vicd)-mu_c*met(akg))/vflux;
    flux(vakgd) = ConvinienceKinetics(model,pvec,met,find(vakgd));
    knidx(vakgd) = 1;
    %vsucoa
    [~,vflux] = ConvinienceKinetics(model,pvec,met,find(vsuca));
    pvec.Vmax(vsuca) = (-flux(vakgd)+mu_c*met(sucoa))/vflux;
    flux(vsuca) = ConvinienceKinetics(model,pvec,met,find(vsuca));
    knidx(vsuca) = 1;
    %vmdh
%     [~,vflux] = ConvinienceKinetics(model,pvec,met,find(vmdh));
%     pvec.Vmax(vmdh) = (-flux(vme1)+flux(vfum)+flux(vmals)-mu_c*met(mal))/vflux;
%     flux(vmdh) = ConvinienceKinetics(model,pvec,met,find(vmdh));
%     knidx(vmdh) = 1;
    %vpyk
    [~,vflux] = ConvinienceKinetics(model,pvec,met,find(vpyk));
    pvec.Vmax(vpyk) = (flux(veno)-flux(vppc)-flux(vpts)-mu_c*met(pep))/vflux;
    flux(vpyk) = ConvinienceKinetics(model,pvec,met,find(vpyk));
    knidx(vpyk) = 1;
    
    %testing to see whether thermodynamics is the problem
%     pvec.Vmax(1:model.nt_rxn) = 1;
    %solve model to obtain concentrations
    [model,batch,solverP] = initModel(model,150);
    variable.MC = met;
    [~,batch] = initConcentration(model,batch,variable);
    [Sol,finalSS,status] =...
    callODEsolver(model,pvec,variable,[],batch,solverP);
    %rxnlist
    
    rxnlist = {''};
    mu_c = biomass_flux(model,met,[],flux);
    
    
end

function [flux,pvec,knidx] = linearfluxes(model,flux,met,mu_c,knidx,pvec,metid,rxnid)
    vmet_f = model.S(metid,:);
    vknw = setdiff(find(vmet_f),find(rxnid));
    [~,vflux] = ConvinienceKinetics(model,pvec,met,find(rxnid));
    pvec.Vmax(rxnid) = (model.S(metid,vknw)*flux(vknw)-mu_c*met(metid))/vflux;
    knidx(rxnid) = 1;
    flux(rxnid) = ConvinienceKinetics(model,pvec,met,find(rxnid));
return

function [flux,pvec,knidx] = pathway_pp(model,S,pvec,flux,met,mu_c,bmrx,knidx)
vpgi = strcmpi(model.rxns,'pgi');
vgnd = strcmpi(model.rxns,'gnd');
vpfk = strcmpi(model.rxns,'pfk');
vfbp = strcmpi(model.rxns,'fbp');
vrpi = strcmpi(model.rxns,'rpi');
vtkt1 = strcmpi(model.rxns,'tkt1');
vtkt2 = strcmpi(model.rxns,'tkt2');
vtal = strcmpi(model.rxns,'tala');

f6p = strcmpi(model.mets,'f6p[c]');
ru5p = strcmpi(model.mets,'ru5p-D[c]');
r5p = strcmpi(model.mets,'r5p[c]');
x5p = strcmpi(model.mets,'xu5p-D[c]');
s7p = strcmpi(model.mets,'s7p[c]');
e4p = strcmpi(model.mets,'e4p[c]');

[vind,~] = find([vpfk vtal vtkt2 vrpi vtkt1]);
[~,vflux] = ConvinienceKinetics(model,pvec,met,vind);

%vpfk   
pvec.Vmax(vpfk) = -(2*flux(vgnd)-flux(vpgi)-mu_c*...
            (-3*S(e4p,1)+2*S(r5p,1)-S(f6p,1)-3*met(e4p)-...
            met(f6p)+2*met(r5p)+2*met(ru5p)-4*met(s7p)+2*met(x5p)))/...
            vflux(vind==find(vpfk));            
%vtala
pvec.Vmax(vtal) = (-flux(vgnd)-mu_c*(-S(r5p,1)+S(e4p,1)+met(e4p)-met(r5p)-...
                                met(ru5p)+2*met(s7p)-met(x5p)))/...
                                vflux(vind==find(vtal));
%vtkt2
pvec.Vmax(vtkt2) = (-flux(vgnd)-mu_c*(-S(r5p,1)+2*S(e4p,1)+2*met(e4p)-...
                                met(r5p)-met(ru5p)+2*met(s7p)-met(x5p)))/...
                                vflux(vind==find(vtkt2));
%vrpi
pvec.Vmax(vrpi) = (-flux(vgnd)+mu_c*(met(ru5p)-S(e4p,1)-met(e4p)-met(s7p)+met(x5p)))/...
              vflux(vind==find(vrpi));
%vrpe
%pvec.Vmax(vrpe) = (mu_c*(met(x5p)-met(e4p)-met(s7p)-S(e4p,1)))/vflux(vind==find(vrpe));          
%vtkt1
pvec.Vmax(vtkt1) = (flux(vgnd)+mu_c*(S(e4p,1)-S(r5p,1)+met(e4p)-met(r5p)-...
                               met(ru5p)+met(s7p)-met(ru5p)))/...
                               vflux(vind==find(vtkt1));   

knidx(vind) = 1;                           
flux(vind) = ConvinienceKinetics(model,pvec,met,vind);

return

%everything else