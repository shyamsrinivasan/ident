% solveEqilibriumODE
% script for running ODEs and MATCONT for multiple sets of parameters
% regardless of the source of parameter samples
% should be run from Kotte2014_Script
[allxdyn,allxeq,allfdyn,allfeq] =...
 solveODEonly(npts,M,model,allpvec,opts,tspan,...
              allxdyn,allxeq,allfdyn,allfeq);

for ipt = 1:npts
    fprintf('\nIteration #%d of %d Equilibrium Continuation...\n',ipt,npts);
    
    % ap = 9; - fix it in the calling function from here on out
    xeq = allxeq(:,ipt);
    pvec = allpvec(ipt,:);
    
    % run MATCONT
    [data,y,p] =...
    execMATCONT(@KotteMATCONT,@Kottecont_fluxcalc,xeq,pvec,ap,fluxg,model);
    if ~isempty(data) && size(data.s1,1)>2
%         bifurcationPlot(data.flux,data.s1,data.f1,[5,3]);
        hbif1 =...
        bifurcationPlot(data.x1,data.s1,data.f1,[4,1],@getKotteaxislabels,ap);    
%         bifurcationPlot(data.x1,data.s1,data.f1,[4,2]);
%         bifurcationPlot([data.flux;data.x1(end,:)],data.s1,data.f1,[6,5]);
    end
    
    % save MATCONT results
    s.(['pt' num2str(ipt)]) = data;
    
    % get the mss for y and p
    if ~isempty(data)
        [yss,iyval,fyval] = parseMATCONTresult(data.s1,y);
        [pss,ipval,fpval] = parseMATCONTresult(data.s1,p);
        [fss,ifval,ffval] = parseMATCONTresult(data.s1,data.flux);
    end
        
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
        nLP = size(s1,1)-2; % [initial lp1 lp2 final] 
        if nLP > 0
            fprintf('Vector %d has %d Limit Points\n',ipt,nLP);
            mssid = union(mssid,ipt);
            nss(ipt) = nLP;
        end
    else
        fprintf('No convergence at %d\n',ipt);
    end
end