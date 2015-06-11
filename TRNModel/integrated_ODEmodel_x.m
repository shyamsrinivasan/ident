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
Mref = data.Yref(ng(1)+ng(2)+ng(3)+2:sum(ng)-ng(5)-ng(6)+1);%intracellular metabolites
CXref = data.Yref(sum(ng)-ng(5)-ng(6)+2:sum(ng)-ng(6)+1);%Complexes
%% %Variables
% mRNA = Y(1:ng(1));%umole
%concentration in mmole
protind = ng(1)+(1:ng(2));
cmind = ng(1)+PMind_R;
mt_ind = ng(1)+ng(2)+ng(3)+2:sum(ng)-ng(5)-ng(6)+1;
% cx_ind = sum(ng)-ng(5)-ng(6)+2:sum(ng)-ng(6)+1;
% RegCMPLX_P = model.RegCMPLX_P;
% RegCMPLX_M = model.RegCMPLX_M;
mx_ind = sum(ng)-ng(6)+2:sum(ng)+1;
regcn = [Y(protind).*Pref;...%protein umoles
         Y(cmind).*PMref;...%common protein umoles
         Y(ng(1)+rnapind).*RNAPref;...%RNAP umoles
         Y(mt_ind).*Mref;...%metabolite mmoles
         %Y(cx_ind).*CXref;...%Complexes ?units?
         data.ext_MC];%extracellular metabolite mmoles

%All Metabolites
Mabs = [regcn(cmind-ng(1)).*1e-3;...%mmole
        regcn([mt_ind,mx_ind]-ng(1))];%mmole
% Mabs = regcn([cmind,mt_ind,mx_ind]-ng(1));%mmole 
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
    %Coefficients to determine binding affinity scaled to umole and mmole 
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
end 
%% %Proteins & Protein Complexes
%Transcription only proteins
%Translation - Proteins and Protein Complexes
% Tlation = beta*(model.ptrate.*model.trate)*(Y(1:ng(1)).*mRNAref)./data.Yref(ng(1)+1:end);%umole/gDCW.s
Tlation = beta*(model.trate)*(Y(1:ng(1)).*mRNAref)./data.Yref(ng(1)+1:end);%umole/gDCW.s
CPLXConc = fun(model,[regcn(protind-ng(1));...%umole
                      regcn(cmind-ng(1));...%umole
                      regcn(rnapind);...%umole
                      regcn(rnapind+1:end).*1e3]);%umole
                  
dYdt(protind) = Tlation(protind-ng(1)) +...%umole/gDCW.s
                CPLXConc(protind-ng(1)) -...%umole/gDCW.s
                (pdecay+gr_flux).*regcn(protind-ng(1));%umole/gDCW.s
            

%dYdt(noRegCLX) = Trasnlation(noRegCLX) + model.S(noReg_M,:)*flux - decay

%RNAP ploymerase - constant concentration
dYdt(ng(1)+rnapind) = 0;    
%% Complexes
% dYdt(cx_ind) = model.Kb(cx_ind-ng(1)).*regcn(RegCMPLX_P).*regcn(RegCMPLX_M) -...
%                model.Kub(cx_ind-ng(1)).*regcn(cx_ind-ng(1));
          
%% %Fluxes 
%All Protein concentrations
EC = Y([protind,cmind,ng(1)+rnapind]).*[Pref;PMref;RNAPref]*1e-3;%from umole to mmole
% EC = Y(Vind+ng(1)).*Pref(Vind)*1e-3;%from umole to mmole
flux(Vind) = ConvinienceKinetics(model,model.pmeter,Mabs,bm_rxn,Vind,EC);
if ~isempty(model.Vupind)
    %flux(model.Vupind) = PTSuptake();
    flux(model.Vupind) = data.vuptake;%mmole/gDCW.s %Vuptake = Vdemand
end
if ~isempty(model.Vexind)
%     flux(model.Vexind) = 0.001*Mabs(intm_ind);
    flux(model.Vexind) = (1000*EC(model.Vexind)).*Mabs(intm_ind);%mmole/gDCW.s
end
if ~isempty(model.Vex)
    %mmole/gDCW.s
    flux(model.Vex) = ConvinienceKinetics(model,model.pmeter,Mabs,bm_rxn,model.Vex,EC);
end
if ~isempty(model.VnoEnz)
    %mmole/gDCW.s    
    flux(model.VnoEnz) = func1(model,model.VnoEnz,Mabs);%noEnzFlux();
end  
% flux(5) = 0;
flux(bm_rxn) = gr_flux;%s-1
%Print Fluxes during every call to function
%printflux(flux,fluxind);
%Biomass
% dYdt(end) = (gr_flux)*Y(end);

%Transcription and metabolic proteins (Eg. PTS)
% dYdt(ng(1)+PMind_R) = model.trate(PMind_R,:)*(Y(1:ng(1)).*mRNAref)./PMref +...
%                       model.S(PMind_M,:)*flux -...
%                       (pdecay + gr_flux).*Y(ng(1)+PMind_R);
%% %Metabolites(Intracelllular)
metindx = ng(3)+1:ng(3)+ng(4);%mt_ind-ng(1)-ng(2)-ng(3)-1;
mt_ind = setdiff(mt_ind-ng(1),model.noRegCLX)+ng(1);
metindx = setdiff(metindx,model.noRegCLX_M);
mindx = setdiff(1:ng(4),model.noRegCLX_M-ng(3));
dYdt(mt_ind) = model.S(metindx,:)*flux./Mref(mindx) -...%mmole/gDCW.s
               Mabs(metindx)*gr_flux;                     %mmole/gDCW.s
           
%Transcription and metabolic proteins (Eg. PTS)
dYdt(cmind) = Tlation(PMind_R) +...%umole/gDCW.s
              CPLXConc(PMind_R) +...%umole/gDCW.s
              model.S(PMind_M,:)*(flux.*1e3) -...%mmole/gDCW.s -> umole/gDCW.s
              (pdecay + gr_flux).*regcn(cmind-ng(1));%umole/gDCW.s           

dYdt(model.noRegCLX+ng(1)) = model.S(model.noRegCLX_M,:)*flux./Mref(model.noRegCLX_M-ng(3)) -...
                             (pdecay+gr_flux).*Mabs(model.noRegCLX_M);                      
           
        
           
%Fixing perturbed variables to zero
if isfield(data,'pbind')
    dYdt(data.pbind) = 0;
end
% dYdt = changeC(t,EC,Mrel,dYd t,model,ng,gr_flux);
% if any(Y(Y<0))
%     flag = 1;
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