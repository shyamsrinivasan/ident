function [allxdyn,allxeq,allfdyn,allfeq,slope,alltout] =...
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
alltout = zeros(length(tspan),npts);
                 

for ipt = 1:npts
    fprintf('\nIteration #%d of %d Equilibrium Integration...',ipt,npts);
    
    % change in pvec
    pvec = allpvec(ipt,:);
    
    % new equilibrium solution
    givenModel = @(t,x)KotteODE(t,x,model,pvec);
    [tout,yout] = ode45(givenModel,tspan,ival,opts);
    
    if size(tout,1) == size(tspan,2)
        allxdyn(:,:,ipt) = yout';
        alltout(:,ipt) = tspan';
    else
        allxdyn(:,1:size(tout,1),ipt) = yout';
        alltout(1:size(tout,1),ipt) = tout;
    end
    allxeq(:,ipt) = yout(end,:)';   
    
    % slope of ODE trajectory (dx/dt) at each tout using ADMAT        
    slope = zeros(size(yout,2),length(tout));
    for it = 1:length(tout)
        slope(:,it) = givenModel(tout(it),yout(it,:)');        
    end
    
    % calculation of fluxes for allxeq
    if size(tout,1)==size(tspan,2)
        allfdyn(:,:,ipt) =...
        Kotte_givenFlux([allxdyn(:,:,ipt);...
                        repmat(model.PM,1,size(allxdyn,2))],pvec,model); 
    else
        allfdyn(:,1:size(tout,1),ipt) =...
        Kotte_givenFlux([allxdyn(:,:,ipt);...
                        repmat(model.PM,1,size(allxdyn,2))],pvec,model); 
    end
    allfeq(:,ipt) = Kotte_givenFlux([allxeq(:,ipt);model.PM],pvec,model); 
    
    
    % optional  - plot information
%      hf = plotKotteVariables(tout,yout,1);
%     plotKotteVariables(tout,allfdyn(:,:,ipt)',2);
    drawnow
%     close(hf);
    fprintf('Complete\n');
end