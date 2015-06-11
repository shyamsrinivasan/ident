function [rr,flag,newdata] = integratedDAEmodel(t,Y,YP,data,model,FBAmodel)
%% %Necessary Quantities
ng = data.ng;
nvar = model.nvar;%%[mRNA,Protein,IntraMetabolite,ExtraMetabolite]'
Vind = model.Vind;
[intm_ind,~] = find(model.S(:,model.Vexind)<0);
rnapind = model.rnap_ind;
CLXid = model.CLX;
CLXA = model.CLX_A;
CLXB = model.CLX_B;
CLXact = model.CLXact;

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
%0 - Algebraic, 1- Differential Variable
varind = data.varind;
var1s = find(varind);
varind(var1s) = 0;
varind(setdiff(1:length(varind),var1s)) = 1;
%Indices in Y or YP
pid = ng(1)+(1:ng(2));
cmid = ng(1)+PMind_R;
mtid = ng(1)+ng(2)+ng(3)+2:sum(ng)-ng(5)-ng(6)+1;
mxid = sum(ng)-ng(6)+2:sum(ng)+1;

%Concentrations
%All Regulators - umole/gDCW
regcn = [Y(pid).*Pref;...               %protein umoles
         Y(cmid).*PMref;...             %common protein umoles
         Y(ng(1)+rnapind).*RNAPref;...  %RNAP umoles
         Y(mtid).*Mref;...              %metabolite mmoles         
         data.ext_MC];                  %extracellular metabolite mmoles

%All Metabolites - mmole/gDCW
Mabs = [Y(cmid).*PMref.*1e-3;...%mmole
        Y(mtid).*Mref;...       %mmole
        Y(mxid)];               %mmole
    
%Indices without cxind
pidx = setdiff(pid,find(varind));
cmidx = setdiff(cmid,find(varind));
mtidx = setdiff(mtid,find(varind));
mxidx = setdiff(mxid,find(varind));

svarind = ismember(model.Metabolites,model.Regulators(find(varind)-ng(1)));
csidx = setdiff(PMind_M,find(svarind));

cxid = varind;
cxidx = find(cxid)-ng(1);

%ODE Residual System
rr = zeros(nvar,1);
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
    RNAP = 3e-5;%umole/gDCW
    rr(igene) = RNAP*(alpha/brate(igene) + alphac*pactnr/(1+pactdr))/mRNAref(igene)-...    
                Y(igene)*(mdecay + gr_flux) - YP(igene); %umole/gDCW.s
end 
%% %Proteins & Protein Complexes
%Transcription only proteins
%Translation - Proteins and Protein Complexes
% Tlation = beta*(model.ptrate.*model.trate)*(Y(1:ng(1)).*mRNAref)./data.Yref(ng(1)+1:end);%umole/gDCW.s
Tlation = beta*(model.trate)*(Y(1:ng(1)).*mRNAref)./data.Yref(ng(1)+1:end);%umole/gDCW.s
%Complexes Regulatory and Nonregulatory and active non active regulator
%forms
[AEres,DEres] = complexAE(model,cxidx,[Y(pid);...          %umole
                                       Y(cmid);...         %umole
                                       Y(rnapind+ng(1));...%umole
                                       Y(mtid).*1e3;...    %umole
                                       Y(mxid).*1e3]);     %umole
%Phosphorylated proteins                                   
[AEres,DEres] = SIGres(model,cxidx,AEres,DEres,[Y(pid).*1e-3;...          %mmole
                                                Y(cmid).*1e-3;...         %mmole
                                                Y(rnapind+ng(1)).*1e-3;...%mmole
                                                Y(mtid);...
                                                Y(mxid)]);                %mmole);                                   

%Residuals for Protein concentrations
rr(pidx) = Tlation(pidx-ng(1)) +...             %umole/gDCW.s              
           DEres(pidx-ng(1)) -...               %umole/gDCW.s
           (pdecay+gr_flux).*Y(pidx)-YP(pidx);  %umole/gDCW.s
       
%Algebraic Variables and corresponding residuals
rr(logical(cxid)) = AEres;                      %umole/gDCW

%RNAP ploymerase - constant concentration
rr(ng(1)+rnapind) = -YP(ng(1)+rnapind);    
         
%% %Fluxes 
%All Protein concentrations
EC = Y([pid,cmid,ng(1)+rnapind]).*[Pref;PMref;RNAPref]*1e-3;%from umole to mmole
% EC = Y(Vind+ng(1)).*Pref(Vind)*1e-3;%from umole to mmole
flux(Vind) = ConvinienceKinetics(model,model.pmeter,Mabs,bm_rxn,Vind,EC);
if ~isempty(model.Vupind)
    %flux(model.Vupind) = PTSuptake();
    flux(model.Vupind) = [0.005 0];%data.vuptake%mmole/gDCW.s %Vuptake = Vdemand
end
if ~isempty(model.Vexind)
%     flux(model.Vexind) = 0.001*Mabs(intm_ind);
    flux(model.Vexind) = (1000*EC(model.Vexind)).*Mabs(intm_ind);%mmole/gDCW.s
end
if ~isempty(model.Vex)
    %mmole/gDCW.s
    flux(model.Vex) = ConvinienceKinetics(model,model.pmeter,Mabs,bm_rxn,model.Vex,EC);
end
% if ~isempty(model.VnoEnz)
% %     mmole/gDCW.s    
%     flux(model.VnoEnz) = func1(model,model.VnoEnz,Mabs);%noEnzFlux();
% end  
% flux(5) = 0;
flux(bm_rxn) = gr_flux;%s-1
%Print Fluxes during every call to function
%printflux(flux,fluxind);
%Biomass
% dYdt(end) = (gr_flux)*Y(end);

%Proteins that act as Metabolites in S
rr(cmidx) = Tlation(cmidx-ng(1))+...                %umole/gDCW.s
            DEres(cmidx-ng(1))+...                  %umole/gDCW.s
            model.S(csidx,:)*(flux.*1e3)-...        %umole/gDCW.s
            (pdecay + gr_flux).*Y(cmidx)-YP(cmidx); %umole/gDCW.s
            

%Transcription and metabolic proteins (Eg. PTS)
% dYdt(ng(1)+PMind_R) = model.trate(PMind_R,:)*(Y(1:ng(1)).*mRNAref)./PMref +...
%                       model.S(PMind_M,:)*flux -...
%                       (pdecay + gr_flux).*Y(ng(1)+PMind_R);
%% %Metabolites(Intracelllular)
%Intracellular metabolites
metind = setdiff(ng(3)+1:ng(3)+ng(4),find(varind)-ng(1)-ng(2)-1);
rr(mtidx) = model.S(metind,:)*flux./Mref(1:length(mtidx)) +...%mmole/gDCW.s
            DEres(mtidx-ng(1)).*1e-3-...                      %mmole/gDCW.s
            Y(mtidx)*gr_flux - YP(mtidx);                     %mmole/gDCW.s
           
%Transcription and metabolic proteins (Eg. PTS)
% rr(cmidx) = Tlation(PMind_R) +...%umole/gDCW.s
%             DEres(PMind_R) +...%umole/gDCW.s
%             model.S(PMind_M,:)*(flux.*1e3) -...%mmole/gDCW.s -> umole/gDCW.s
%             (pdecay + gr_flux).*regcn(cmidx-ng(1))-YP(cmidx);%umole/gDCW.s           

% rr(model.noRegCLX+ng(1)) = model.S(model.noRegCLX_M,:)*flux./Mref(model.noRegCLX_M-ng(3)) -...
%                              (pdecay+gr_flux).*Mabs(model.noRegCLX_M)-...
%                              YP(model.noRegCLX+ng(1));                      
           
        
           
%Fixing perturbed variables to zero
if isfield(data,'pbind')
    rr(data.pbind) = 0;
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
  
return