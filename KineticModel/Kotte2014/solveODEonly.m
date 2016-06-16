function [allxdyn,allxeq,allfdyn,allfeq,slope] =...
         solveODEonly(npts,ival,model,allpvec,opts,tspan,...
                      allxdyn,allxeq,allfdyn,allfeq)

if nargin<10     
    allfeq =...
    zeros(length(Kotte_givenFlux([ival;model.PM],allpvec(1,:),model)),npts);
end
if nargin<9
    allfdyn =...
    zeros(length(Kotte_givenFlux([ival;model.PM],allpvec(1,:),model)),...
          length(tspan),npts);
end
if nargin<8
    allxeq = zeros(length(ival),npts);
end
if nargin<7
    if length(tspan)>2
        allxdyn = zeros(length(ival),length(tspan),npts);
    end        
end
                 

for ipt = 1:npts
    fprintf('\nIteration #%d Equilibrium Integration...',ipt);
    
    % change in pvec
    pvec = allpvec(ipt,:);
    
    % new equilibrium solution
    givenModel = @(t,x)KotteODE(t,x,model,pvec);
    [tout,yout] = ode45(givenModel,tspan,ival,opts);
    if ~exist('allxdyn','var')
        allxdyn = zeros(length(ival),length(tout),npts);
    end
    allxdyn(:,:,ipt) = yout';
    allxeq(:,ipt) = yout(end,:)';   
    
    % slope of ODE trajectory (dx/dt) at each tout using ADMAT        
    slope = zeros(size(yout,2),length(tout));
    for it = 1:length(tout)
        slope(:,it) = givenModel(tout(it),yout(it,:)');        
    end
%     slope = [];
    
%     Jxact = KottegivenJacobian(ival,pvec,model);
    
%     dX = zeros(length(tout),
%     dX = givenModel(
    
%     xeq = yout(end,:)';
    
    % calculation of fluxes for allxeq
    allfeq(:,ipt) = Kotte_givenFlux([allxeq(:,ipt);model.PM],pvec,model); 
    for itout = 1:length(tout)
        allfdyn(:,itout,ipt) = Kotte_givenFlux([allxdyn(:,itout,ipt);model.PM],pvec,model); 
    end
    
    % optional  - plot information
     plotKotteVariables(tout,yout,1);
%     plotKotteVariables(tout,allfdyn(:,:,ipt)',2);

    fprintf('Complete\n');
end