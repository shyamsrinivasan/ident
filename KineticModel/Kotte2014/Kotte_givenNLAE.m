function dM = Kotte_givenNLAE(kmrgd,model,pvec)

d = pvec(13);
dM = zeros(3,1);
if ~isempty(model)
    allmc = [kmrgd;model.PM];
else
    allmc = kmrgd;
end
flux = Kotte_givenFlux(allmc,pvec,model);
% differential equations
% PEP
dM(1) = flux(1) - flux(4);
% FBP
dM(2) = flux(4) - flux(3);
% enzymes
% E
dM(3) = flux(2) - d*allmc(3);