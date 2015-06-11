function [dYdt,flag,newdata] = trnODE_ADMAT(Y,Extra)
model = Extra.trnmodel;
FBAmodel = Extra.FBAmodel;
data = Extra.data;

%Necessary Quantities
ng = data.ng;
nreg = sum(ng(2:end));
nvar = sum(ng)+1;%there are ng(2) fluxes
SSval = data.SSval;

%Variables
mRNA = Y(1:ng(1));
regulator = Y(ng(1)+1:ng(1)+nreg);

%Parameter Initialization
reg_full = SSval(ng(1)+1:end-1).*(1+regulator);
metab_full = SSval(ng(1)+ng(2)+1:sum(ng)).*(1+regulator(ng(2)+1:sum(ng)-ng(1)));
% prot_full = [regulator(1:ng(2)).*model.maxReg(1:ng(2));regulator(ng(2)+1:end)];
% prot_full= regulator.*model.maxReg;
[Coefficient,brate] = parameter_return(data.par,model);
alpha = 0.01;%transcription rate

% Pss %ss protein concentration
pdecay = data.pdecay;
mdecay = data.mdecay;
% V_ind = [FBAmodel.Vin_ind,FBAmodel.Vout_ind];
% Vint = FBAmodel.Vind;
%%
%ODE System
dYdt = zeros(nvar,1);
%mRNA Transcription
for igene = 1:ng(1)
    %Binding Affinity = Protein Concentration/Binding Constant
    regind = logical(model.RS(igene,:));
    bindaff = reg_full(regind)./Coefficient(igene,regind)';
    prot_name = model.Regulators(regind);    
%     bindaff = (allProt(regind).*model.maxProt(regind))./Coefficient(igene,regind)';
    [pactnr,pactdr] =  regactivity(prot_name,...
                       model.RS(igene,regind),...
                       model.GeneRules{igene},...
                       bindaff,...
                       data);
    %freg(igene) = pactnr/(1+pactdr);
    %dG/dt = rst*(pactnr/(1+pactdr)) - decay*G;
    dYdt(igene) = alpha*SSval(end)/SSval(igene)*(1+Y(end))*...
                  (1/brate(igene) + pactnr/(1+pactdr))- mdecay*(1+mRNA(igene));    
    %rst = dRNAP/dt*(1/brate + pactnr/(1+pactdr));
    %dYdt(ngene+igene) = dYdt(2*ngene+nprot+1)*(1 + pactnr/(1+pactdr));
end 
%Proteins
dYdt(ng(1)+1:ng(1)+ng(2)) = (SSval(1:ng(1))./SSval(ng(1)+1:ng(1)+ng(2))).*(model.trate(1:ng(2),:)...
                            *Y(1:ng(1)))-(1+Y(ng(1)+1:ng(1)+ng(2)))*pdecay;
%Flux = [Protein]
flux = 0.01*SSval(ng(1)+1:ng(1)+ng(2)).*(1+Y(ng(1)+1:ng(1)+ng(2)));
% Vmax = Y(ng(1)+1:ng(1)+ng(2)).*model.maxReg(1:ng(2));
% MC = Y(ng(1)+ng(2)+1:ng(1)+ng(2)+ng(3));
% [mind,~] = find(FBAmodel.S(1:ng(3),Vint)<0);
% flux(Vint) = Vmax(Vint).*(MC(mind)./([5;5;5]+MC(mind)));
% Y(ng(1)+ng(2)+ng(3)+ng(4)+1:end-1) = flux;
%Metabolites(Intracelllular)
dYdt(ng(1)+ng(2)+1:ng(1)+ng(2)+ng(3)) = (1./SSval(ng(1)+ng(2)+1:ng(1)+ng(2)+ng(3))).*...
                                        ((FBAmodel.S(1:ng(3),:)*flux) -...
                                        0.2*SSval(end)*(1+Y(end)).*...
                                        (1+regulator(ng(2)+1:ng(2)+ng(3))));
%Metabolites(Extracellualr)
% vsupply = flux(V_ind);
dYdt(ng(1)+ng(2)+ng(3)+1:ng(1)+ng(2)+ng(3)+ng(4)) =...
                                (1./SSval(ng(1)+ng(2)+ng(3)+1:sum(ng))).*...
                                ([0;0]-...
                                FBAmodel.M*FBAmodel.S(1:ng(3),:)*flux)-...
                                (1+regulator(ng(2)+ng(3)+1:end));
%External ODE = Supply - FBAmodel.M*FBAmodel.S*flux;
%Fluxes
% dflux = dYdt(ng(1)+1:ng(1)+ng(2)).*model.maxReg(1:ng(2));
% dflux = zeros(ng(2),1);
% Dr = [5;5;5]+MC(mind); 
% dVmax = dYdt(ng(1)+1:ng(1)+ng(2)).*model.maxReg(1:ng(2));
% dMetab = dYdt(ng(1)+ng(2)+1:ng(1)+ng(2)+ng(3));
% term1 = (MC(mind)./Dr).*dVmax(Vint);
% term2 = (Vmax(Vint)./Dr).*dMetab(mind);
% dflux(Vint) = term1 + term2 - Vmax(Vint).*MC(mind).*dYdt(mind)./(Dr.^2);
% dflux(~Vint) = dVmax(~Vint);
% dYdt(sum(ng)+1:sum(ng)+ng(2)) = dflux;
%RNAP
dYdt(end) = (1/SSval(end))*sum((metab_full(ng(3)+1:ng(3)+ng(4)).*...
            [1e-5;0])./([50;1]+metab_full(ng(3)+1:ng(3)+ng(4))));
%Variables = mRNA, rst, regProteins, recProteins, eRNAP, Metab, Metabxt,
%Biomass


    
%Regulators = [Proteins;Metabolites]
%Metab = Y(
% Metab = Y(2*ngene+nprot+2:2*ngene+nprot+1+nmetab);
% Metabxt = Y(2*ngene+nprot+nmetab+2:2*ngene+nprot+1+2*nmetab);
% mCell = Y(end);
% eRNAP = Y(2*ngene+nprot+1);
% allProt = Y(2*ngene+1:2*ngene+nprot);
% %totalProt = sum(Y(2*ngene+1:2*ngene+nprot));
% regProt = Y(2*ngene+1:2*ngene+nregp);
% recProt = Y(2*ngene+nregp+1:2*ngene+nregp+nrecp);

%RNAP
% dYdt(2*ngene+nprot+1) = sum(model.rnapsyn.*Metab./(model.Ksmetab + Metab)) -...
%                         eRNAP*(pdecay + dYdt(end)/mCell);
% freg = zeros(ngene,1);  

%Extracellular Metabolites
%dMxtidt = -(kcat*model.transmat*allProt).*Mxti i = 1 to nmetab
% dYdt(2*ngene+nprot+nmetab+2:2*ngene+nprot+1+2*nmetab) =...
%      -(.08*model.transmat*allProt).*Metabxt;

%Intracellular Metabolites
%dMidt = uptake - yieldi*growth i = 1 to nmetab
% dYdt(2*ngene+nprot+2:2*ngene+nprot+1+nmetab) =...
%     ((0.08*model.transmat*allProt).*Metabxt - model.myield*dYdt(end))./model.molwt;

%Biomass mCell - Based on Intracellular Glucose Concentration
%dMc/dt = [mumax1;mumax2]'*[S1/Ks1+S1;S2/Ks2+S2]
% dYdt(end) = mCell*sum(model.mumax.*Metab./(model.Ksmetab + Metab));

%Protein
%d[regProt;recProt]/dt = [trate, metabRS]*[mRNA;Metabolites] -
%Protein*(decay + growth)
% dYdt(2*ngene+1:2*ngene+nregp) = model.trate(1:nregp)*Y(1:ngene) -...
%                                 Y(2*ngene+1:2*ngene+nregp)*...
%                                 (pdecay + dYdt(end-1)/mCell);
% dYdt(2*ngene+nregp+1:2*ngene+nregp+nrecp) = model.metabRS(nregp+1:nprot)*...
%      (model.Kmax.*Metabxt./(model.Ks + Metabxt)) -...
%      Y(2*ngene+nregp+1:2*ngene+nprot)*(pdecay + dYdt(end-1)/mCell);              


flag = 0;
newdata = [];    
end