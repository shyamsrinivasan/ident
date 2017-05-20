function dM = noisyNLAE_kotte(x,pstruct)
% 
if isfield(pstruct,'p')
    p = pstruct.p;
end
if isfield(pstruct,'model')
    model = pstruct.model;
else
    model = [];
end

% x = x + rand(3,1);
d = p(10);
dM = zeros(3,size(x,2));
if ~isempty(model)
%     PM = cons(model.PM,x);
    PM = model.PM;
    allmc = [x;repmat(PM,1,size(x,2))];
else
    allmc = x;
end
% dM = cons(dM,allmc);
flux = noisyflux_kotte(allmc,pstruct);
% differential equations
% PEP
dM(1,:) = flux(1,:) - flux(4,:) - flux(5,:);
% FBP
dM(2,:) = flux(4,:) - flux(3,:);
% enzymes
% E
dM(3,:) = flux(2,:) - d*allmc(3,:);

% generate noise
% noise = rand(3,1);
% dM = dM + repmat(noise,1,size(dM,2));

