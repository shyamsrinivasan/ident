function [dYdt,flag,newdata] = matbalance(t,Y,data,model,FBAmodel)
%Necessary Quantities
ng = data.ng;
nreg = sum(ng(2:end));
nvar = sum(ng)+2;%[mRNA,Protein,Metabolite,Biomass,RNAP]'
nt_rxn = model.nt_rxn;
SSval = data.SSval;
Vind = model.Vind;
Vic_exind = model.Vic_exind;
Vxc_exind = model.Vxc_exind;
bm_ind = model.bm_ind;
Mbio = model.S(:,bm_ind);
Mbio_ind = find(Mbio<0);
%Variables
mRNA = Y(1:ng(1));
regulator = Y(ng(1)+1:ng(1)+nreg);
metab = regulator(ng(2)+1:sum(ng)-ng(1));
%Parameters

% reg_full = SSval(ng(1)+1:sum(ng)).*(1+regulator);
% metab_full = SSval(ng(1)+ng(2)+1:sum(ng)).*(1+regulator(ng(2)+1:sum(ng)-ng(1)));
[Coefficient,brate] = parameter_return(data.par,model);
alpha = 1e-4;%s-1 transcription rate
kcat = 100;%s-1
% D = 1;%Dilution rate h-1
pdecay = data.pdecay;
mdecay = data.mdecay;

%ODE System
dYdt = zeros(nvar,1);
% dYdt = cons(dYdt,Y);
flux = zeros(nt_rxn,1);
%mRNA Transcription
for igene = 1:ng(1)
    %Binding Affinity = Protein Concentration/Binding Constant
    regind = logical(model.RS(igene,:));
    bindaff = regulator(regind)./Coefficient(igene,regind)';
%     bindaff = reg_full(regind)./Coefficient(igene,regind)';
    prot_name = model.Regulators(regind);    
    [pactnr,pactdr] =  regactivity(prot_name,...
                       model.RS(igene,regind),...
                       model.GeneRules{igene},...
                       bindaff,...
                       data); 
    dYdt(igene) = alpha*Y(end)*(1/brate(igene)+pactnr/(1+pactdr))-mdecay*mRNA(igene);
%     dYdt(igene) = alpha*SSval(end)/SSval(igene)*(1+Y(end))*...
%                   (1/brate(igene) + pactnr/(1+pactdr))- mdecay*(1+mRNA(igene));        
end 
%Fluxes 
EC = Y(Vind+ng(1))*1e-3;
flux(Vind) = ConvinienceKinetics(model,model.pmeter,metab(1:ng(3)),bm_ind,EC);
% flux(Vind) = kcat*Y(Vind+ng(1))*1e-3;
% flux(Vind) = kcat*SSval(Vind+ng(1)).*(1+Y(Vind+ng(1)));
flux(Vic_exind) = kcat*Y(Vic_exind+ng(1))*1e-3;%Exchange(Intracellular)
% flux(Vic_exind) = kcat*SSval(Vic_exind+ng(1)).*(1+Y(Vic_exind+ng(1)));
flux(Vxc_exind) = -flux(Vic_exind);%Exchange(Extracellular)

gr_flux = sum([0.8;0.6].*metab(Mbio_ind)./([.8;.1]+metab(Mbio_ind)))/3600;
flux(bm_ind) = gr_flux*Y(end-1);
D = gr_flux*3600;
% flux(bm_ind) = gr_flux*Y(end-1)/3600;%Vbiomass
%flux(bm_ind) = -sum(model.S(Mbio_ind,bm_ind).*dYdt(Mbio_ind+ng(1)+ng(2)));
% flux(bm_ind) = -sum(model.S(Mbio_ind,bm_ind).*SSval(Mbio_ind+ng(1)+ng(2)).*dYdt(Mbio_ind+ng(1)+ng(2)));
% gr_flux = flux(bm_ind)*Y(end-1);
%gr_flux = gr_max*prod(metab_full(Mbio_ind)./([.8;.1]+metab_full(Mbio_ind)));
% flux(bm_ind) = gr_flux*(1+Y(end-1))*SSval(end-1)/(1e-3*bm_molwt);
%Proteins
dYdt(ng(1)+1:ng(1)+ng(2)) = model.trate(1:ng(2),:)*Y(1:ng(1))-(pdecay+gr_flux)*Y(ng(1)+1:ng(1)+ng(2));
% dYdt(ng(1)+1:ng(1)+ng(2)) = (SSval(1:ng(1))./SSval(ng(1)+1:ng(1)+ng(2))).*...
%                             (model.trate(1:ng(2),:)*(1+Y(1:ng(1))))-...
%                             (pdecay+gr_flux)*(1+Y(ng(1)+1:ng(1)+ng(2)));
%Metabolites(Intracelllular)
dYdt(ng(1)+ng(2)+1:ng(1)+ng(2)+ng(3)) = model.S(1:ng(3),:)*flux-gr_flux*metab(1:ng(3));
% dYdt(ng(1)+ng(2)+1:ng(1)+ng(2)+ng(3)) =...
% (1./SSval(ng(1)+ng(2)+1:ng(1)+ng(2)+ng(3))).*(model.S(1:ng(3),:)*flux);%-...
% (gr_flux).*(1+regulator(ng(2)+1:ng(2)+ng(3)));
%Biomass
dYdt(end-1) = (gr_flux - D/3600)*Y(end-1);
% dYdt(end-1) = (gr_flux - D)*(1+Y(end-1));
%Metabolites(Extracellular) %No Equation for Extracellular
% dYdt(ng(1)+ng(2)+ng(3)+1:sum(ng)) = D/3600*([0;0]-regulator(ng(2)+ng(3)+1:end))-...
%                                     FBAmodel.S(ng(3)+1:ng(3)+ng(4),:)*flux;
%                                     FBAmodel.M*FBAmodel.S(1:ng(3),:)*flux;
%Metabolites(Extracellular)
dYdt(ng(1)+ng(2)+ng(3)+1:sum(ng)) = D/3600*([10;10]-regulator(ng(2)+ng(3)+1:end))-...
                                    FBAmodel.M*FBAmodel.S(1:ng(3),:)*flux;
% dYdt(ng(1)+ng(2)+ng(3)+1:sum(ng)) =...
% (1./SSval(ng(1)+ng(2)+ng(3)+1:sum(ng))).*(D*[0;0]-FBAmodel.M*FBAmodel.S(1:ng(3),:)*flux)-...
% D.*(1+regulator(ng(2)+ng(3)+1:end));
%RNA Polymerase
% dYdt(end) = (9e-3)*(55.4e-2)*Y(end-1)-(pdecay+gr_flux)*Y(end);
dYdt(end) = .554*gr_flux*0.005*(1/450e-3)-(pdecay+gr_flux)*Y(end);
% dYdt(end) = ((9e-3)*(55.4e-2)*SSval(end-1)/(SSval(end)*450))*Y(end-1) -...
%             (pdecay+gr_flux)*(1+Y(end));

flag = 0;
newdata = [];                
  
end