function [mc,assignFlag,delGr,vCorrectFlag] = assignConc(x,model,bounds)
%assign concentrations to the model at large
if size(x,2)>1
    mc = zeros(length(model.mets),size(x,2));
else
    mc = zeros(length(model.mets),1);
end

%appending mets with same lb and ub
x = [x;repmat(bounds.x,1,size(x,2))];
bounds.A = [bounds.A(:,1:length(bounds.mets)) bounds.A_kn];
bounds.mets = [bounds.mets;bounds.mets_kn];

%check for delGr values
delGr = getdelGr(bounds,x);

vCorrectFlag = zeros(length(bounds.rxns),size(x,2));
vCorrectFlag(sign(delGr.*repmat(bounds.Vss,1,size(x,2)))<0) = 1;

assignFlag = zeros(length(model.mets),1);
for im = 1:length(bounds.mets)
    tfm = strcmpi(bounds.mets{im},model.mets);
    if any(tfm)
        mc(tfm,:) = x(im,:);
        assignFlag(tfm) = 1;
    end
end

