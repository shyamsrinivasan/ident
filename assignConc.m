function [mc,assignFlag] = assignConc(x,model,bounds)
%assign concentrations to the model at large
if size(x,2)>1
    mc = zeros(length(model.mets),size(x,2));
else
    mc = zeros(length(model.mets),1);
end
assignFlag = zeros(length(model.mets),1);
for im = 1:length(bounds.mets)
    tfm = strcmpi(bounds.mets{im},model.mets);
    if any(tfm)
        mc(tfm,:) = x(im,:);
        assignFlag(tfm) = 1;
    end
end