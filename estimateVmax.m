function estimateVmax(model,pvec,met)
bmrx = model.bmrxn;
Vuptake = model.Vuptake;
S = -model.S(:,bmrx);
%known fluxes
knflux = {'G6PDH2r','pfk','pps','ppc','pdh','me1'};
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
vfba = strcmpi(model.rxns,'fba');
vtpi = strcmpi(model.rxns,'tpi');
%metabolites
g6p = strcmpi(model.mets,'g6p[c]');
pgl = strcmpi(model.mets,'6pgl[c]');
pgc = strcmpi(model.mets,'6pgc[c]');
f6p = strcmpi(model.mets,'f6p[c]');
fdp = strcmpi(model.mets,'fdp[c]');
dhap = strcmpi(model.mets,'dhap[c]');
g3p = strcmpi(model.mets,'g3p[c]');
dpg = strcmpi(model.mets,'13dpg[c]');
%unknown fluxes
unkflux = {'glcpts','pgi','pgl','gnd','fbp','fba','rpe','rpi','tkt1',...
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
    %fbp,rpe,rpi,tkt1,tala,tkt2
    [flux,knidx] = pathway_pp(S,flux,met,knidx);
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

function [flux,knidx] = pathway_pp(S,flux,met,knidx)
vgnd = strcmpi(model.rxns,'gnd');
vpfk = strcmpi(model.rxns,'pfk');
vfbp = strcmpi(model.rxns,'fbp');
vrpe = strcmpi(model.rxns,'rpe');
vrpi = strcmpi(model.rxns,'rpi');
vtkt1 = strcmpi(model.rxns,'tkt1');
vtkt2 = strcmpi(model.rxns,'tkt2');
vtal = strcmpi(model.rxns,'tala');

ru5p = strcmpi(model.mets,'ru5p-D[c]');
r5p = strcmpi(model.mets,'r5p[c]');
x5p = strcmpi(model.mets,'xu5p-D[c]');
s7p = strcmpi(model.mets,'s7p[c]');
e4p = strcmpi(model.mets,'e4p[c]');

%vfbp    
flux(vfbp) = 2*flux(vgnd)+flux(vpfk)-v(pgi)-mu_c*...
            (-3*S(e4p,bmrx)+2*S(r5p,bmrx)-S(f6p,bmrx)-3*met(e4p)-...
            met(f6p)+2*met(r5p)+2*met(ru5p)-4*met(s7p)+2*met(x5p));
            
%vtala
flux(vtal) = -flux(vgnd)-mu_c*(-S(r5p,1)+S(e4p,1)+met(e4p)-met(r5p)-...
                                met(ru5p)+2*met(s7p)-met(xu5p));
%vtkt2
flux(vtkt2) = -flux(vgnd)-mu_c*(-S(r5p,1)+2*S(e4p,1)+2*met(e4p)-...
                                met(r5p)-met(ru5p)+2*met(s7p)-met(xu5p));
%vrpe
flux(vrpe) = mu_c*(met(xu5p)-met(e4p)-met(s7p)-S(e4p,1));
%vrpi
flux(vrpi) = -flux(vgnd)+mu_c*(met(ru5p)-S(e4p,1)-met(e4p)-met(s7p)+met(x5p));
%vtkt1
flux(vtkt1) = flux(vgnd)+mu_c*(S(e4p,1)-S(r5p,1)+met(e4p)-met(r5p)-...
                               met(ru5p)+met(s7p)-met(ru5p));   

knidx([vfbp vtal vtkt2 vrpe vrpi vtkt1]) = 1;                           
flux([vfbp vtal vtkt2 vrpe vrpi vtkt1]) =...
ConvinienceKinetics(model,pvec,met,find([vfbp vtal vtkt2 vrpe vrpi vtkt1]));

return

%everything else