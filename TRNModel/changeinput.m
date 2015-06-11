function [MetabConc] = changeinput(model,metabname,metabconc,MetabConc)

if nargin < 4
    MetabConc = zeros(15,1);
end

if nargin < 3
    fprintf('Value >= 0 Required for Input Metabolite Concentration');
    return
end

if length(metabname) ~= length(metabconc)
    fprintf('Dimensions of Metabolite & Concentrations do not match');
    return
end

if length(metabname) <= 1
    metabindx = strcmp(metabname,model.MetabName);
    MetabConc(metabindx) = metabconc;
else
    for imetab = 1:length(metabname)
        metabindx = strcmp(metabname{imetab},model.MetabName);
        MetabConc(metabindx) = metabconc(imetab);
    end
end

end