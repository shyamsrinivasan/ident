function mc = separate_slack(x,model,bounds)
if ~isempty(x)
    mc = x(1:length(bounds.mets));
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


