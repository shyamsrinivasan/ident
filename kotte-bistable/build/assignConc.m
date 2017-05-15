function [mc,assignFlag,delGr,vCorrectFlag] = assignConc(x,model,bounds)
%assign concentrations to the model at large
if size(x,2)>1
    mc = zeros(length(model.mets),size(x,2));
else
    mc = zeros(length(model.mets),1);
end

%appending mets with same lb and ub
x = [x;repmat(bounds.x_kn,1,size(x,2))];
bounds.A = [bounds.A bounds.A_kn];
bounds.mets = [bounds.mets;bounds.mets_kn];

%check for delGr values
delGr = getdelGr(bounds,x);

vCorrectFlag = zeros(length(bounds.rxns),size(x,2));
vCorrectFlag(sign(delGr.*repmat(bounds.Vss,1,size(x,2)))<=0) = 1;

assignFlag = zeros(length(model.mets),1);
for im = 1:length(bounds.mets)
    tfm = strcmpi(bounds.mets{im},model.mets);
    if any(tfm)
        mc(tfm,:) = x(im,:);
        assignFlag(tfm) = 1;
    end
end
assignFlag = logical(assignFlag);
vCorrectFlag = logical(vCorrectFlag);

if ~isempty(delGr)
    delGr = assignRxns(delGr,model,bounds);
end

%assign same concentrations to intra and extracellular co2,o2,h,pi and h2o
%if available
% met = {'co2','o2','h','pi','h2o'};
% for imet = 1:length(met)
%     tfbnd = strcmpi(bounds.mets,[met{imet} '[c]']);
%     tfmod = strcmpi(model.mets, [met{imet} '[e]']);
%     if any(tfbnd) && any(tfmod)
%         mc(tfmod,:) = x(tfbnd,:);
%         assignFlag(tfmod) = 1;
%     end    
% end

