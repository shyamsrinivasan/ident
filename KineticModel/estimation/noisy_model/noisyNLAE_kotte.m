function dM = noisyNLAE_kotte(x,p)

d = p(10);
dM = zeros(4,size(x,2));

% generate noisy flux
flux = kotte_flux_noCAS(x,p);
flux = flux + 2*rand(5,1);

% differential equations
% PEP
dM(1,:) = flux(1,:) - flux(4,:) - flux(5,:);
% FBP
dM(2,:) = flux(4,:) - flux(3,:);
% enzymes
% E
dM(3,:) = flux(2,:) - d*x(3,:);
% acetate
dM(4,:) = x(4,:) - x(4,:);


% generate noise
% noise = rand(3,1);
% dM = dM + repmat(noise,1,size(dM,2));

