function [AEres,DEres] = SIGres(model,cxind,AEres,DEres,regcn)
VnoEnz = model.VnoEnz;
AEflux = zeros(length(VnoEnz),1);
for ine = 1:length(VnoEnz)
    subsind = model.S(:,VnoEnz(ine)) < 0;   
    prodind = model.S(:,VnoEnz(ine)) > 0;
    [subrid,prodrid] = getindex(model,subsind,prodind);
    %PTS Phosphoryl transfer reactions
    if any(strcmpi(model.Regulators(subrid),'EI')) &&...
        any(strcmpi(model.Regulators(subrid),'pep[c]'))
        %EI + pep <---> pyr + P~EI
        %Parameters 
        kI = 1e17;%gDCW/mmole.s
        kuI = 6.67e16;%gDCW/mmole.s
        AEflux(ine) = kI*prod(regcn(subrid)) -...
                      kuI*prod(regcn(prodrid));               %mmole/gDCW.s
    elseif any(strcmpi(model.Metabolites(subsind),'P_EI')) &&...
        any(strcmpi(model.Metabolites(subsind),'HPr'))   
        %P~EI + HPr <---> EI + P~HPr
        %Parameters 
        kII = 3.33e17;%gDCW/mmole.s
        kuII = 1.33e17;%gDCW/mmole.s
        AEflux(ine) = kII*prod(regcn(subrid)) -...
                      kuII*prod(regcn(prodrid));              %mmole/gDCW.s
    elseif any(strcmpi(model.Metabolites(subsind),'P_HPr')) &&...
        any(strcmpi(model.Metabolites(subsind),'EIIA'))
     %P~Hpr + EIIA <---> HPr + P~EIIA 
     %Parameters 
     kIII = 6.1e13;%gDCW/mmole.s from Rowler et al., 2000
     kuIII = 4.7e13;%gDCW/mmole.s from Rowler et al., 2000
     AEflux(ine) = kIII*prod(regcn(subrid)) -...
                   kuIII*prod(regcn(prodrid));                %mmole/gDCW.s
    end        
end

res = model.S(:,VnoEnz)*(AEflux.*1e3);                        %umole/gDCW.s
p_ind = ~cellfun('isempty',regexp(model.Regulators,'^(?:P_)\w'));
px_ind = intersect(find(p_ind),cxind);
cmid = ismember(model.Metabolites,model.Regulators(px_ind));

AEres(ismember(cxind,px_ind)) = res(cmid);                    %umole/gDCW
deind = setdiff(1:length(DEres),cxind);
for kid = 1:length(deind)
    tfr = strcmpi(model.Regulators{deind(kid)},model.Metabolites);
    if any(tfr)
        DEres(deind(kid)) = res(tfr);                         %umole/gDCW.s
    end
end
return
function [subrid,prodrid] = getindex(model,subsind,prodind)
subs = model.Metabolites(subsind);
prud = model.Metabolites(prodind);
subrid = false(length(model.Regulators),1);
prodrid = false(length(model.Regulators),1);
for is = 1:length(find(subsind))
    tfr = strcmpi(subs{is},model.Regulators);
    if any(tfr)
        subrid(tfr) = 1;
    end        
end
for ip = 1:length(find(prodind))
    tfr = strcmpi(prud{ip},model.Regulators);
    if any(tfr)
        prodrid(tfr) = 1;
    end        
end
return