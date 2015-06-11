function [F] = admat_intODEmodel(Y,defparval,ng,model)
%% %Necessary Quantities
nvar = sum(ng)-ng(4);%%[mRNA,Protein,IntraMetabolite,ExtraMetabolite]'
% nt_rxn = model.nt_rxn;
Vind = model.Vind;
Vic_exind = model.Vic_exind;
Vxc_exind = model.Vxc_exind;
[~,Vuptake] = find(model.S(:,Vic_exind)>0);
[int_ind,Vexcrt] = find(model.S(:,Vic_exind)<0);
Vuptake = Vic_exind(Vuptake);
Vexcrt = Vic_exind(Vexcrt);
bm_ind = model.bm_ind;
%% %Parameters
[Coefficient,brate] = parameter_return(model.allpar,model);
gmax = model.gmax;%h-1
pdecay = defparval.pdecay;
mdecay = defparval.mdecay;
model.Yref = ones(nvar,1);
mRNAref = model.Yref(1:ng(1));
Pref = model.Yref(ng(1)+1:ng(1)+ng(2));
Mref = model.Yref(ng(1)+ng(2)+1:sum(ng)-ng(4));
%% %Variables
% Y = abs(Y);
% mRNA = Y(1:ng(1));%umole
%[protein;int_metab;ext_metab];
regulator =...
[Y(ng(1)+1:ng(1)+ng(2)).*Pref*1e-3;...
 Y(ng(1)+ng(2)+1:sum(ng)-ng(4)).*Mref;...
 model.ext_MC];%mmole
Mabs = regulator(ng(2)+1:sum(ng)-ng(1));%mmole
Mrel = Mabs(1:ng(3))./Mref;%mmole
%ODE System
F = zeros(nvar,1);
F = cons(F,Y);
% flux = zeros(FBAmodel.nt_rxn,1);
flux = model.Vss;
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
    bindaff = regulator(regind)./Coefficient(igene,regind)';
    prot_name = model.Regulators(regind);    
    [pactnr,pactdr] =  regactivity(prot_name,...
                       model.RS(igene,regind),...
                       model.GeneRules{igene},...
                       bindaff,...
                       defparval); 
    RNAP = 3e-5;%uM/gDCW
    F(igene) = RNAP*(alpha/brate(igene) + alphac*pactnr/(1+pactdr))/mRNAref(igene)-...    
                  Y(igene)*(mdecay + gr_flux);
%     dYdt(igene) = alpha*RNAP*(1/brate(igene)+pactnr/(1+pactdr))/mRNAref(igene)-...
%                   Y(igene)*(mdecay + gr_flux);    
end 
%% %Proteins
F(ng(1)+1:ng(1)+ng(2)) = (0.06*beta*model.trate(1:ng(2),:))*...
                            (Y(1:ng(1)).*mRNAref./Pref)-...
                            (pdecay+gr_flux).*Y(ng(1)+1:ng(1)+ng(2));
%% %Fluxes 
EC = Y(Vind+ng(1)).*Pref(Vind)*1e-3;
flux(Vind) = ConvinienceKinetics(model,model.pmeter,Mabs(1:ng(3)),bm_ind,EC);
flux(Vuptake) = 5.56e-1;%mmole/gDCW.s data.Vuptake;%Vuptake = Vdemand
% flux(Vuptake) = 10*(Y(Vuptake+ng(1))*1e-3.*Pref(Vuptake))./...
%                 (1e-2+Y(Vuptake+ng(1))*1e-3.*Pref(Vuptake)).*Mabs(ext_ind);
% flux(Vuptake) = (1000000*Y(Vuptake+ng(1))*1e-3.*Pref(Vuptake)).*...
%                 ((Mabs(ext_ind)./(Mabs(ext_ind)+K)));%-...
% %                 (Mabs(min_log)./(Mabs(min_log)+[100])));
% flux(Vexcrt) = (10000*Y(Vexcrt+ng(1))*1e-3.*Pref(Vexcrt)).*...
%                (Mabs(int_ind)./(100+Mabs(int_ind)));
% flux(Vexcrt) = -1e-4;
flux(Vexcrt) = (10000*Y(Vexcrt+ng(1))*1e-3.*Pref(Vexcrt)).*Mabs(int_ind);
% flux(Vexcrt) = (kcat*Y(Vexcrt+ng(1))*1e-3.*Pref(Vexcrt)).*...
%                (Mabs(int_ind)./(100+Mabs(int_ind)));%-...
%                Mabs(mex_log)./(100+Mabs(mex_log)));
flux(Vxc_exind) = flux(Vic_exind);%Exchange(Extracellular)
flux(bm_ind) = gr_flux;
% flux(Vexcrt) = model.S(int_ind,~V_nexcrt)*flux(~V_nexcrt);
%Biomass
% dYdt(end) = (gr_flux)*Y(end);
%% %Metabolites(Intracelllular)
F(ng(1)+ng(2)+1:ng(1)+ng(2)+ng(3)) = model.S(1:ng(3),:)*flux./Mref(1:ng(3))-...
                                        Mrel(1:ng(3))*gr_flux;  
end