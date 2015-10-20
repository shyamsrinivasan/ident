function delGr = getdelGr(model,mc)
A = model.A(:,1:length(model.mets));
if size(mc,2)>1
    delGr = zeros(length(model.rxns),size(mc,2));
else
    delGr = zeros(length(model.rxns),1);
end

for ipt = 1:size(mc,2)
    Gamma = A*mc(:,ipt);
    for i = 1:length(model.rxns)
        delGr(i,ipt) = 0.008314*298*(Gamma(i)-log(model.Keq(i)));
        if abs(delGr(i,ipt))<1e-8
            delGr(i,ipt) = 0;
        end
    end
end