function dM = noisyNLAE_kotte(x,model,p)

d = p(10);
dM = zeros(3,size(x,2));
if ~isempty(model)
%     PM = cons(model.PM,x);
    allmc = [x;repmat(PM,1,size(x,2))];
else
    allmc = x;
end
dM = cons(dM,allmc);
flux = noisyflux_kotte(allmc,p,model);
% differential equations
% PEP
dM(1,:) = flux(1,:) - flux(4,:) - flux(5,:);
% FBP
dM(2,:) = flux(4,:) - flux(3,:);
% enzymes
% E
dM(3,:) = flux(2,:) - d*allmc(3,:);