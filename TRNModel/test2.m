function [Solution] = test2(trnmodel,ngene,tnreg,ngap,InConc,t0,initval)

if nargin < 7 || length(initval) ~= ngene+tnreg
    initval = ones(ngene+tnreg,1);
    initval = initval/(1E-15*6.023E+23);
end

AbsTol = zeros(size(initval));
AbsTol(AbsTol==0) = 1.e-6;
options = CVodeSetOptions('LMM','Adams',...
                          'NonlinearSolver','Functional',...
                          'RelTol',1.e-8,...
                          'AbsTol',AbsTol);
CVodeInit(@matbalance,'Adams','Functional',t0,initval,options);

%for i = 1:6000
    [status,t,dY] = CVode(t0,'Normal');
    stats = CVodeGetStats;
    fprintf('t = %0.2e  order = %1d step = %0.2e',t,stats.qlast,stats.hlast);
    if status == 0
        t0 = t0*i;
    end
%end
stats = CVodeGetStats;
CVodeFree;

                          

%[t,dY] = ode23(@matbalance,tstart,initval);
Solution.time = t;
Solution.variable = dY;
Solution.g = dY(:,1:ngene);
Solution.p = dY(:,ngene+1:end);


function dYdt = matbalance(t,Y)
decay = 0.00016; %Decay Rate in s-1
G = Y(1:ngene);
P = Y(ngene+1:ngene+ngap);
Pm = Y(ngene+ngap+1:end);

dYdt = zeros(ngene+tnreg,1);
tP = [P;Pm];

for igene = 1:ngene
    %binding affinity is determined on the basis of the same concentrations
    %i.e the initial concentrations. Change in protein concentration as a
    %result of the solution to the ODE is not reflected in subsequent steps
    %of the ODE solution.This has been corrected w/ tP
    bindaff = bindaffinity(trnmodel,trnmodel.RS(igene,:),trnmodel.GeneRules{igene},igene,tP);
    pact = singlepromoteractivity(trnmodel,trnmodel.RS(igene,:),trnmodel.GeneRules{igene},bindaff,igene);    
    dYdt(igene) = pact - decay*G(igene);
end

% for iprot = ngene+1:ngene+ngap
%     dYdt(iprot) = trnmodel.trate(iprot,:)*G - decay*P(iprot);
dYdt(ngene+1:ngene+ngap) = trnmodel.trate*G - decay*P;
% end
prodrate = recpprod(InConc);
dYdt(ngene+ngap+1:end) = prodrate - decay*Pm;

%Material balance function
% dGdt(ngene x 1) = pact(ngene x 1) - decay*G(ngene x 1);
% dPdt(ngap x 1) = translationrate(ngap x ngene)*G(ngene x 1) - decay*P(ngap x 1);
% dPmdt(tnreg-ngap x 1) = rate(tnreg-ngap x nmetab)*M(nmetab x 1) -
% decay*Pm(tnreg-ngap x 1);

end
function [prod] = recpprod(InConc)
%prodrate - Function similar to pactivity.m that delivers a tnreg-ngap x 1 
%matrix of production rates for metabolite TCS-based TFs. This rate
%function will be based on standard reversible MM kinetics
    nmetab = length(trnmodel.Metabolite);
    %# Receptors = # of metabolites (Signalling molecules correspond to different metabolites)      
    if nmetab == tnreg-ngap
        prod = zeros(tnreg-ngap,1);
        for iregp = 1:nmetab
            prod(iregp) = trnmodel.Kmax(ngap+iregp,iregp)*InConc(iregp)/...
                (trnmodel.Ks(ngap+iregp,iregp)+InConc(iregp));
        end
    end   
    
end
end