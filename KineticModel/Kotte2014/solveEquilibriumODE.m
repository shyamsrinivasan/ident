% solveEqilibriumODE
% script for running ODEs and MATCONT for multiple sets of parameters
% regardless of the source of parameter samples
% should be run from Kotte2014_Script
[allxdyn,allxeq,allfdyn,allfeq] =...
 solveODEonly(npts,M,model,allpvec,opts,tspan,...
              allxdyn,allxeq,allfdyn,allfeq);

for ipt = 1:npts
    fprintf('\nIteration #%d Equilibrium Continuation...\n',ipt);
    
    ap = 9;
    xeq = allxeq(:,ipt);
    pvec = allpvec(ipt,:);
    
    % run MATCONT
    [data,y,p] = execMATCONT(xeq,pvec,ap,fluxg,model);
    bifurcationPlot(data.flux,data.s1,data.f1,[5,3]);
    bifurcationPlot(data.x1,data.s1,data.f1,[4,1]);
    bifurcationPlot([data.flux;data.x1(end,:)],data.s1,data.f1,[6,5]);
    
    % save MATCONT results
    s.(['pt' num2str(ipt)]) = data;
    
    % get the mss for y and p
    [yss,iyval,fyval] = parseMATCONTresult(data.s1,y);
    [pss,ipval,fpval] = parseMATCONTresult(data.s1,p);
    [fss,ifval,ffval] = parseMATCONTresult(data.s1,data.flux);
        
    fprintf('Equilibrium Continuation Complete\n');
end

% using fsolve to determine equilibrium points    
%     gfun = @(x)Kotte_givenNLAE(x,model,pvec);
%     options = optimoptions('fsolve','Display','iter',...
%                                     'TolFun',1e-12,'TolX',1e-12,...
%                                     'MaxFunEvals',10000,...
%                                     'MaxIter',5000);
%     [xf,fval,exitflag,output,jacobian] = fsolve(gfun,M,options);


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