function estimateVmax(model,pvec,met)
bmrx = model.bmrxn;
Vuptake = model.Vuptake;
S = -model.S(:,bmrx);
%known fluxes
knflux = {'G6PDH2r','pdh','me1'};
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
%metabolites
g6p = strcmpi(model.mets,'g6p[c]');
pgl = strcmpi(model.mets,'6pgl[c]');
pgc = strcmpi(model.mets,'6pgc[c]');
ru5p = strcmpi(model.mets,'ru5p-D[c]');
fdp = strcmpi(model.mets,'fdp[c]');
dhap = strcmpi(model.mets,'dhap[c]');
g3p = strcmpi(model.mets,'g3p[c]');
dpg = strcmpi(model.mets,'13dpg[c]');
%unknown fluxes
unkflux = {'glcpts','pgi','pgl','gnd','fba','rpe','rpi','tkt1',...
           'tala','tkt2','tpi','gapd','pgk','pgm','eno','ACONTa','cs',...
           'fum','sucdi','mals','icl','icdh','akgdh','SUCOAS',...
           'mdh','ppc','pyk'};  
%initial mu guess
mu_g = 0.1;
%initial mu calculated
mu_c = biomass_flux(model,met,[],Vuptake);
flux(bmrx) = mu_c;
Vmax = zeros(model.nt_rxn,1);

while abs(mu_c-mu_g)>=1e-4
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
    
    
    %rxnlist
    
    rxnlist = {''};
    
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