% solveEqilibriumODE
% script for running ODEs and MATCONT for multiple sets of parameters
% regardless of the source of parameter samples
% should be run from Kotte2014_Script
for ipt = 1:npts
    fprintf('Iteration #%d Equilibrium Integration...',ipt);
    % change in pvec
    pvec = allpvec(ipt,:);
%     pvec(idp) = allpvec(ipt,:);
    
    % new equilibrium solution
    givenModel = @(t,x)KotteODE(t,x,model,pvec);
    [tout,yout] = ode45(givenModel,tspan,M,opts);
    allyoutss(:,ipt) = yout(end,:)';    
%     plotKotteVariables(tout,yout,1);    
    
    gfun = @(x)Kotte_givenNLAE(x,model,pvec);
    options = optimoptions('fsolve','Display','iter',...
                                    'TolFun',1e-12,'TolX',1e-12,...
                                    'MaxFunEvals',10000,...
                                    'MaxIter',5000);
%     [xf,fval,exitflag,output,jacobian] = fsolve(gfun,M,options);
    xeq = allyoutss(:,ipt);
    allxeq(:,ipt) = xeq;
%     allxf(:,ipt) = xf;
%     allflag(1,ipt) = exitflag;
%     xeq = xf;
    
    fprintf('Complete\n');
    
    % continuation from initial equilibrium - initialization
    fprintf('Iteration #%d Equilibrium Continuation...\n',ipt);
    % run MATCONT
    runMATCONT
    
    % save MATCONT results
    s.(['pt' num2str(ipt)]).s1 = s1;
    s.(['pt' num2str(ipt)]).x1 = x1;
    s.(['pt' num2str(ipt)]).f1 = f1;
    
    fprintf('Equilibrium Continuation Complete\n');
end

% check which solutions have mss
mssid = [];
for ipt = 1:npts
    s1 = s.(['pt' num2str(ipt)]).s1;
    nLP = size(s1,1);
    if nLP > 2
        fprintf('Vector %d has %d Steady States\n',ipt,nLP);
        mssid = union(mssid,ipt);
    end
end