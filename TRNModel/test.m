function [Solution] = test(trnmodel,prodrate,ngene,tnreg,ngap)
clc
initval = ones(ngene+tnreg,1);
initval = initval/(1E-15*6.023E+23);


[t,dGdt] = ode23(@matbalance,[0,100],initval);


function dGdt = matbalance(t,G)
decay = 0.00016; %Decay Rate in s-1
%translation = 0.0006; %Translation rate in s-1
%ProtConc = trnmodel.ProtConc;

protconc = G(ngene+1:ngene+ngap);
[bindaff] = getbindingaffinity(trnmodel,ngene,protconc);
[pact] = pactivity(trnmodel,bindaff,ngene);


%Material balance function
% dGdt(ngene x 1) = pact(ngene x 1) - decay*G(ngene x 1);
% dPdt(ngap x 1) = translationrate(ngap x ngene)*G(ngene x 1) - decay*P(ngap x 1);
% dPmdt(tnreg-ngap x 1) = rate(tnreg-ngap x nmetab)*M(nmetab x 1) -
% decay*Pm(tnreg-ngap x 1);

dGdt = zeros(ngene+tnreg,1);
dGdt(1:ngene) = pact - decay*G(1:ngene);
dGdt(ngene+1:ngene+ngap) = trnmodel.trate*G(1:ngene) - decay*G(ngene+1:ngene+ngap);
dGdt(ngene+gap+1:end) = prodrate - decay*G(ngene+gap+1:end);

%prodrate - Function similar to pactivity.m that delivers a tnreg-ngap x 1 
%matrix of production rates for metabolite TCS-based TFs. This rate
%function will be based on standard reversible MM kinetics
end


end
