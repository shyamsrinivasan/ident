function dM = Kotte_CNLAE(kmrgd,model,pvec,ssflux)

if nargin<4
    d = pvec(13);
    dM = zeros(3,1);
    allmc = [kmrgd;model.PM];
    % substitute with Convenience Kinetics
    flux = Kotte_CFlux(allmc,pvec,model);    
    % differential equations
    % PEP
    dM(1) = flux(1) - flux(4);
    % FBP
    dM(2) = flux(4) - flux(3);
    % enzymes
    % E
    dM(3) = flux(2) - d*allmc(3);
else
    allmc = [kmrgd;model.PM];
    % substitute with Convenience Kinetics
    flux = Kotte_CFlux(allmc,pvec,model); 
    dM = zeros(4,1);
    dM(1) = ssflux(1) - flux(1);
    dM(2) = ssflux(2) - flux(2);
    dM(3) = ssflux(3) - flux(3);
    dM(4) = ssflux(4) - flux(4);
end