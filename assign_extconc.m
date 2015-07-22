% function [mConc] = assign_extconc(name,conc,model)
% Function to assign pre-determined constant external
% metabolite concentrations
function [mConc] = assign_extconc(name,conc,model)
    mConc = zeros(model.next_metab,1);
    exter_mind = ~cellfun('isempty',regexp(model.mets,'\w(?:\[e\])$'));
    exter_metab = model.mets(exter_mind);
    if ~isempty(exter_metab)
        if length(name) == length(conc)
            for iname = 1:length(name)
                tfm = strcmpi(name{iname},exter_metab);
                mConc(tfm,1) = conc(iname);
            end
        else
            fprintf('\n# Metabolites & # Concentrations DO NOT Match');
        end
    end
end