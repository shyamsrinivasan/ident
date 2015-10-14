function [mc,assignFlag] = separate_slack(x,model,bounds)
if ~isempty(x)
    mc = exp(x(1:length(bounds.mets)));
end

if ~isempty(mc)
    [mc,assignFlag] = assignConc(mc,model,bounds);        
end


