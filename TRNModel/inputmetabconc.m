%Function to decide input metabolite concentration in the absence of
%absolute concentrations
%function InConc = inputconc(metabs,conc)
function InConc = inputmetabconc(metabName,metabConc,model)
%nmetab = length(trnmodel.Metabolite);
InConc = zeros(length(model.Metabolite),1);
if length(metabName) == length(metabConc)
    for imetab = 1:length(metabName)
        InConc(strcmpi(metabName{imetab},model.Metabolite)) = ...
            metabConc(imetab);
        
        if ~any(strcmpi(metabName{imetab},model.Metabolite))
            fprintf('Specified Metabolite(s) is NOT present in the model\n');
            fprintf('(%d)%s - Not Present\n',imetab,metabName{imetab});
        end        
            
    end
end
end
        
        
