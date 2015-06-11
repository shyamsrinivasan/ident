function [dYdt,flag,newdata] = integrated_ODEmodel(t,Y,data,model,FBAmodel)
%% %Necessary Quantities
ng = data.ng;
nvar = model.nvar;%%[mRNA,Protein,IntraMetabolite,ExtraMetabolite]'
Vind = model.Vind;
[intm_ind,~] = find(model.S(:,model.Vexind)<0);
rnapind = model.rnap_ind;
PMind_R = model.PMind_R;
PMind_M = model.PMind_M;
bm_rxn = model.bmrxn;

Mbio = model.S(:,bm_rxn);
Mbio_ind = find(Mbio<0);
%% %Parameters
[Coefficient,brate] = parameter_return(data.par,model);
% alpha = 1e-4;
% kcat = 10000;%s-1
kcat = data.kcat;%s-1
K = 20;%umole
gmax = data.gmax;%h-1
% D = 0;%Dilution rate h-1
pdecay = data.pdecay;
mdecay = data.mdecay;
mRNAref = data.Yref(1:ng(1));
Pref = data.Yref(ng(1)+1:ng(1)+ng(2));
PMref = data.Yref(ng(1)+ng(2)+1:ng(1)+ng(2)+ng(3));
RNAPref = data.Yref(ng(1)+rnapind);
Mref = data.Yref(ng(1)+ng(2)+ng(3)+2:sum(ng)+1-ng(5));%intracellular metabolites
%% %Variables
% mRNA = Y(1:ng(1));%umole
%concentration in mmole
protind = ng(1)+(1:ng(2));
cmind = ng(1)+PMind_R;
mt_ind = ng(1)+ng(2)+ng(3)+2:sum(ng)-ng(5)+1;
mx_ind = sum(ng)-ng(5)+2:sum(ng)+1;
regcn = [Y(protind).*Pref;...%protein umoles
         Y(cmind).*PMref;...%common protein umoles
         Y(ng(1)+rnapind).*RNAPref;...%RNAP umoles
         Y(mt_ind).*Mref;...%metabolite mmoles
         data.ext_MC];%extracellular metabolite mmoles

% Mabs = [regcn(ng(2):ng(2)+ng(3)-1)
Mabs = regcn([cmind,mt_ind,mx_ind]-ng(1));%mmole - All Metabolites
Mrel = Mabs(ng(3)+1:ng(3)+ng(4))./Mref;%mmole - Intracellular 
%ODE System
dYdt = zeros(nvar,1);
% dYdt = cons(dYdt,Y);
% flux = zeros(FBAmodel.nt_rxn,1);
flux = data.flux;
%Biomass
gr_flux = gmax/3600;%sum(gmax.*Mabs(Mbio_ind)./([50]+Mabs(Mbio_ind)))/3600;%s-1
alphac = 0.06;%s-1
alpha = 40/(233/(gr_flux*3600)^2+78);%s-1 transcription rate
beta = 60/(82.5/(gr_flux*3600)+148);%s-1 translation rate
model.trate(model.trate~=0) = 1;
%% %mRNA Transcription
for igene = 1:ng(1)
    %Binding Affinity = Protein Concentration/Binding Constant
    regind = logical(model.RS(igene,:));
    bindaff = regcn(regind)./Coefficient(igene,regind)';
    prot_name = model.Regulators(regind);    
    [pactnr,pactdr] =  regactivity(prot_name,...
                       model.RS(igene,regind),...
                       model.GeneRules{igene},...
                       bindaff,...
                       data); 
    RNAP = 3e-5;%uM/gDCW
    dYdt(igene) = RNAP*(alpha/brate(igene) + alphac*pactnr/(1+pactdr))/mRNAref(igene)-...    
                  Y(igene)*(mdecay + gr_flux);
%     dYdt(igene) = alpha*RNAP*(1/brate(igene)+pactnr/(1+pactdr))/mRNAref(igene)-...
%                   Y(igene)*(mdecay + gr_flux);    
end 
%% %Proteins
%Transcription only proteins
% protindx = ng(1)+setdiff(1:ng(2),rnapind);
dYdt(protind) = beta*(0.06*model.trate(protind-ng(1),:))*...%(model.ptrate(protind-ng(1),:).*
                (Y(1:ng(1)).*mRNAref)./Pref(protind-ng(1)) -...
                (pdecay+gr_flux).*Y(protind);
%RNAP ploymerase - constant concentration
dYdt(ng(1)+rnapind) = 0;   
                        
%% %Fluxes 
%All Protein concentrations
EC = Y([protind,cmind,ng(1)+rnapind]).*[Pref;PMref;RNAPref]*1e-3;%from umole to mmole
% EC = Y(Vind+ng(1)).*Pref(Vind)*1e-3;%from umole to mmole
flux(Vind) = ConvinienceKinetics(model,model.pmeter,Mabs,bm_rxn,Vind,EC);
if ~isempty(model.Vupind)
    flux(model.Vupind) = data.vuptake;%mmole/gDCW.s %Vuptake = Vdemand
end
if ~isempty(model.Vexind)
%     flux(model.Vexind) = 0.001*Mabs(intm_ind);
    flux(model.Vexind) = (1000*EC(model.Vexind)).*Mabs(intm_ind);
end
if ~isempty(model.Vex)
    flux(model.Vex) = ConvinienceKinetics(model,model.pmeter,Mabs,bm_rxn,model.Vex,EC);
end
% if ~isempty(Vout)
%     flux(Vout) = ConvinienceKinetics(model,model.pmeter,Mabs,bm_rxn,Vout,EC);
% end
% flux(Vuptake) = 5.56e-3;%mmole/gDCW.s data.Vuptake;%Vuptake = Vdemand
% flux(Vuptake) = 10*(Y(Vuptake+ng(1))*1e-3.*Pref(Vuptake))./...
%                 (1e-2+Y(Vuptake+ng(1))*1e-3.*Pref(Vuptake)).*Mabs(ext_ind);
% flux(Vuptake) = (1000000*Y(Vuptake+ng(1))*1e-3.*Pref(Vuptake)).*...
%                 ((Mabs(ext_ind)./(Mabs(ext_ind)+K)));%-...
% %                 (Mabs(min_log)./(Mabs(min_log)+[100])));
% flux(Vexcrt) = (10000*Y(Vexcrt+ng(1))*1e-3.*Pref(Vexcrt)).*...
%                (Mabs(int_ind)./(100+Mabs(int_ind)));
% flux(Vexcrt) = -1e-4;
% flux(Vexcrt) = (100*Y(Vexcrt+ng(1))*1e-3.*Pref(Vexcrt)).*Mabs(int_ind);
% flux(Vexcrt) = (kcat*Y(Vexcrt+ng(1))*1e-3.*Pref(Vexcrt)).*...
%                (Mabs(int_ind)./(100+Mabs(int_ind)));%-...
%                Mabs(mex_log)./(100+Mabs(mex_log)));
% flux(Vxc_exind) = flux(Vic_exind);%Exchange(Extracellular)
flux(bm_rxn) = gr_flux;
% flux(Vexcrt) = model.S(int_ind,~V_nexcrt)*flux(~V_nexcrt);
%Biomass
% dYdt(end) = (gr_flux)*Y(end);

%Transcription and metabolic proteins (Eg. PTS)
dYdt(ng(1)+PMind_R) = model.trate(PMind_R,:)*(Y(1:ng(1)).*mRNAref)./PMref +...
                      model.S(PMind_M,:)*flux -...
                      (pdecay + gr_flux).*Y(ng(1)+PMind_R);
%% %Metabolites(Intracelllular)
metindx = mt_ind-ng(1)-ng(2)-ng(3)-1;
dYdt(mt_ind) = model.S(metindx,:)*flux./Mref(metindx) -...
               Mrel(metindx)*gr_flux;
           
%Fixing perturbed variables to zero
if isfield(data,'pbind')
    dYdt(data.pbind) = 0;
end
% dYdt = changeC(t,EC,Mrel,dYd t,model,ng,gr_flux);
% if any(Y(Y<0))
%     flag = -1;
%     fprintf('Model does not converge\n');
%     newdata = [];
%     return
% else
%     flag = 0;
% end
flag = 0;
newdata = data;
newdata.MC = Mabs(1:ng(3));
newdata.flux = flux;
  
end