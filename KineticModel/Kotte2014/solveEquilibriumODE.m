% solveEqilibriumODE
% script for running ODEs and MATCONT for multiple sets of parameters
% regardless of the source of parameter samples
% should be run from Kotte2014_Script

for ipt = 1:npts
    fprintf('\nIteration #%d Equilibrium Integration...',ipt);
    % change in pvec
    pvec = allpvec(ipt,:);
%     pvec(idp) = allpvec(ipt,:);
    
    if ipt == 1001
        dbstop at 34
    end
    % new equilibrium solution
    givenModel = @(t,x)KotteODE(t,x,model,pvec);
    [tout,yout] = ode45(givenModel,tspan,M,opts);
    allxdyn(:,:,ipt) = yout';
    allxeq(:,ipt) = yout(end,:)';       
    
    gfun = @(x)Kotte_givenNLAE(x,model,pvec);
    options = optimoptions('fsolve','Display','iter',...
                                    'TolFun',1e-12,'TolX',1e-12,...
                                    'MaxFunEvals',10000,...
                                    'MaxIter',5000);
%     [xf,fval,exitflag,output,jacobian] = fsolve(gfun,M,options);
    xeq = yout(end,:)';
     
%     allxf(:,ipt) = xf;
%     allflag(1,ipt) = exitflag;
%     xeq = xf;
    
    % calculation of fluxes for allxeq
    allfeq(:,ipt) = Kotte_givenFlux([allxeq(:,ipt);model.PM],pvec,model); 
    for itout = 1:length(tout)
        allfdyn(:,itout,ipt) = Kotte_givenFlux([allxdyn(:,itout,ipt);model.PM],pvec,model); 
    end
    plotKotteVariables(tout,yout,1);
%     plotKotteVariables(tout,allfdyn(:,:,ipt)',2);
    fprintf('Complete\n');
    
    % continuation from initial equilibrium - initialization
    fprintf('Iteration #%d Equilibrium Continuation...\n',ipt);
    % run MATCONT
%     clear sys x0 v0 ap opt x1 v1 s1 h1 f1
    ap = 9; % index for parameter to be continued on     
%     runMATCONT
    [data,y,p] = execMATCONT(xeq,pvec,ap,fluxg,model);
    
    % save MATCONT results
    s.(['pt' num2str(ipt)]) = data;
    
    % get the mss for y and p
    [yss,iyval,fyval] = parseMATCONTresult(data.s1,y);
    [pss,ipval,fpval] = parseMATCONTresult(data.s1,p);
    [fss,ifval,ffval] = parseMATCONTresult(data.s1,data.flux);
        
    fprintf('Equilibrium Continuation Complete\n');
end

% check which solutions have mss
mssid = [];
nss = zeros(npts,1);
for ipt = 1:npts
    if ~isempty(s.(['pt' num2str(ipt)]))
        s1 = s.(['pt' num2str(ipt)]).s1;
        nLP = size(s1,1);
        if nLP > 2
            fprintf('Vector %d has %d Steady States\n',ipt,nLP);
            mssid = union(mssid,ipt);
            nss(ipt) = nLP;
        end
    else
        fprintf('No convergence at %d\n',ipt);
    end
end