function [vf,assignFlag] = assignRxns(x,model,bounds)
%assign delG and fluxes to the model at large
if size(x,2)>1
    vf = zeros(length(model.rxns),size(x,2));
else
    vf = zeros(length(model.rxns),1);
end

assignFlag = zeros(length(model.rxns),1);
for iv = 1:length(bounds.rxns)
    tfv = strcmpi(bounds.rxns{iv},model.rxns);
    if any(tfv)
        vf(tfv,:) = x(iv,:);        
        assignFlag(tfv) = 1;
    end
end