function [AEres,DEres] = complexAE(model,cxind,regcn)
%Function to calculate Binding and unBinding flux in the formation of
%regulatory complexes like cAMP-CRP
ncx = length(cxind);
ncplx = size(model.CSA,2);
% cxpid = unique(find(model.CSA(:,1:ncplx)));
AEres = zeros(ncx,1);
AssID = zeros(ncx,1);
CSAres = zeros(ncplx,1);
DEres = zeros(size(model.CSA,1),1);
for icp = 1:ncplx
    subind = model.CSA(:,icp) < 0;    
    prodind = model.CSA(:,icp) > 0;
    if any(find(prodind)==cxind)
        if any(strcmpi(model.Regulators{prodind},'cAMP-CRP'))
            %cAMP-CRP -  Malecki et al., 2000
            kc = 9.7;%s-1
            kuc = 0.31;%s-1
            Kanti = 27.5e-15;%umole/gDCW
            %0 = kc/Kanti*[cAMP]*[CRP] - (kuc)*[cAMP-CRP]
            AEres(find(prodind)==cxind) =...
            kc/Kanti*prod(regcn(subind))-kuc*(regcn(prodind));  %umole/gDCW
            CSAres(icp) = AEres(find(prodind)==cxind);          
            AssID(find(prodind)==cxind) = 1;
        end        
        if any(strcmpi(model.Regulators{prodind},'LacI-lac[c]'))
            %LacI
            Klac = 1;%uM Semsey et al., 2013
            tLacI = strcmpi('LacI',model.Regulators);
            tlac = strcmpi('lac[c]',model.Regulators);        
            %0 = LacI(1-1/(1+lac[c]/Klac)) - [LacI-lac[c]]
            AEres(find(prodind)==cxind) =...
            regcn(tLacI)*(1-1/(1+regcn(tlac)/Klac)) - regcn(prodind);%umole/gDCW
            CSAres(icp) = AEres(find(prodind)==cxind);
            AssID(find(prodind)==cxind) = 1;
        end
    end
end

unAssID = cxind(~AssID);
for icp = 1:length(unAssID)
    if any(strcmpi(model.Regulators{unAssID(icp)},'LacI*'))
        %LacI - 
        Klac = 1;%uM Semsey et al., 2013
        tLacI = strcmpi('LacI',model.Regulators);
        tlac = strcmpi('lac[c]',model.Regulators);
        %0 = LacI/(1+lac[c]/Klac) - [LacI*]
        AEres(unAssID(icp)==cxind) =...
        regcn(tLacI)/(1+regcn(tlac)/Klac) - regcn(prodind);     %umole/gDCW
%         CSAres(icp) = AEres(find(prodind)==cxind);
    end    
end


for icp = 1:ncplx
    subind = model.CSA(:,icp) < 0;    
    prodind = model.CSA(:,icp) > 0;
    if any(find(prodind)==cxind)
        if any(strcmpi(model.Regulators{prodind},'cAMP-CRP'))
            DEres(subind) = model.CSA(subind,:)*CSAres;         %umole/gDCW.s
        end
    end
end


% [~,ncrxn] = size(model.CS);
% FormFlux = zeros(ncrxn,1);
% unFormFlux = zeros(ncrxn,1);

%Parameters

% KbcAMP = 2.44e-12;%umole/gDCW - Kremling et al., 2001 Part 2
% % dCMPLXdt = zeros(nreg,1);
% for icrxn = 1:ncrxn
%     subsind = model.CS(:,icrxn)<0;
%     prodind = model.CS(:,icrxn)>0; 
%     tfcAMP = strcmpi(model.Regulators,'cAMP');
%     tfCRP = strcmpi(model.Regulators,'CRP');
%     %if CRP and cAMP
%     if any(strcmpi(model.Regulators(subsind),'CRP')) &&...
%             any(strcmpi(model.Regulators(subsind),'cAMP'))                    
%         kcAMP = kuc + kc*((regcn(tfcAMP))^2/...
%                 (Kanti^2 + 2*Kanti*regcn(tfcAMP) + regcn(tfcAMP)^2));%s-1        
% %         rcAMP = kcAMP*regcn(tfCRP)*(regcn(tfcAMP)/(KbcAMP + regcn(tfcAMP)));%umole/gDCW
%         rcAMP = kcAMP*regcn(tfCRP)*regcn(tfcAMP);%umole/gDCW
%         FormFlux(icrxn) = rcAMP;
%     %Other Complexes
%     else
%         FormFlux(icrxn) = model.Kb(prodind).*prod(regcn(subsind));%umole/gDCW.s
%         unFormFlux(icrxn) = model.Kub(prodind).*prod(regcn(prodind));%umole/gDCW.s
%     end   
% end
% dCMPLXdt = model.CS*(FormFlux - unFormFlux);%umole/gDCW.s
return