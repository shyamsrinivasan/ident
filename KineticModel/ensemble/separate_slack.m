function mc = separate_slack(x,bounds)
if ~isempty(x)
    if size(x,2)>1
        mc = x(:,1:length(bounds.mets));
    else
    mc = x(1:length(bounds.mets));
    end
%     bounds.A = bounds.A(:,1:length(bounds.mets));
%     bounds.lb = bounds.lb(1:length(bounds.mets));
%     bounds.ub = bounds.ub(1:length(bounds.mets));
%     
%     %add known variables to the list
%     mc = [mc;bounds.x];
%     bounds.mets = [bounds.mets;bounds.mets_kn];
%     bounds = rmfield(bounds,{'mets_kn','x'});
end
% assignFlag = [];
% if ~isempty(mc)
%     [mc,assignFlag] = assignConc(mc,model,bounds);        
% end


