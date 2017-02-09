function [allxeq,allfeq] = solveAEonly(npts,ival,model,allpvec,opts,...
                                       allxeq,allfeq)
                  
if nargin<7     
    allfeq =...
    zeros(length(Kotte_givenFlux([ival;model.PM],allpvec(1,:),model)),npts);
end
if nargin<6
    allxeq = zeros(length(ival),npts);
end  

for ipt = 1:npts
    fprintf('\nIteration #%d of %d Nonlinear Algebraic Solution...',ipt,npts);
    
    % change in pvec
    pvec = allpvec(ipt,:);
    
    % new equilibrium solution
    givenModel = @(x)Kotte_givenNLAE(x,model,pvec);
    [xeq,fval,exitflag] = fsolve(givenModel,ival,opts);    
    
    if exitflag
        allxeq(:,ipt) = xeq;   
    end    
    
    % calculation of fluxes for allxeq
    allfeq(:,ipt) = Kotte_givenFlux([allxeq(:,ipt);model.PM],pvec,model);         
    
    fprintf('Complete\n');
end