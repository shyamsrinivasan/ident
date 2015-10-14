function [mc,assignFlag] = separate_slack(x,model,bounds)
if ~isempty(x)
    mc = x(1:length(bounds.mets));
end
assignFlag = [];
% if ~isempty(mc)
%     [mc,assignFlag] = assignConc(mc,model,bounds);        
% end


