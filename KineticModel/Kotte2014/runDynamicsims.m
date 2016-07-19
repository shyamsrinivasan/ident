% runDynamicsims
% simulate from one or more limit points or within their boundaries 
% as initial conditions to test stability/dynamic behaviour of points
% load the relevant simulation dataset
% change to simulate from any given/multiple saddle nodes (points on the
% bistable line)
% load('C:\Users\shyam\Documents\Courses\CHE1125Project\Results\KotteModel\VmaxVariation_July19.mat');
% Plot trajectories and trialing sepratrix vector plots

% change order of parameters to suit new order
% apold = ap;
% alliidpvecold = alliidpvec;
% ap = changeparamsorder('old2newid',apold);
% alliidpvec = changeparamsorder('old2new',alliidpvecold);

% needed variables: alliidpvec,alliidxeq,alliidfeq,tout,ap;
npts = size(alliidpvec,1);
nvar = size(alliidxeq,1);
ndp = size(alliidpvec,3);
tout = tspan;

% FIGmssdynamicswpvec(alliidxdyn,tout,alliidpvec,1,1,1:1000,'conc',...
%                     1:find(tout==100))

opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
colorSpec = chooseColors(4,{'Green','Purple','Red','Orange'});

ac = find(strcmpi(model.mets,'ac[e]'));

% choose nxsspts random points for poseeigind within bounds in
% xLPval
% preliminary evaluation to pre allocate memory for collecting ode sim data
nxsspts = 3;
nsols = zeros(ndp,1);
for idp = 1:ndp
    % changed parameters    
    if isfield(allnss,sprintf('iid%d',idp));
        msspts = find(allnss.(['iid' num2str(idp)]));
        sslps = allnss.(['iid' num2str(idp)])(msspts);
        ss = unique(sslps);
        nss = length(ss);
        % number of different mss
        for iss = 1:nss
            allmsspts = msspts(sslps==ss(iss));
            nmsspts = length(allmsspts);
            % collect corresponding limit points
            for ipt = 1:nmsspts
                index =...
                cat(1,siid.(['iid' num2str(idp)]).(['pt' num2str(allmsspts(ipt))]).s1.index);
                nindexpts = length(index);
                nsols(idp) = nsols(idp)+(nindexpts-1)*nxsspts;
            end
        end
    end
end



% determine the #parameters/combinations that have been changed
isol = 1;
for idp = 1:ndp
    % determine actual # parameters that have been changed
    npar = length(alliidpvec(1,:,idp));
    if npts>1
        diffpar = find(alliidpvec(1,:,idp)~=alliidpvec(2,:,idp));
    end    
    hfig = [];
    hsfig = [];
    hline = [];    
    allival = zeros(nvar,nsols(idp));
    allxeq = zeros(nvar,nsols(idp));
    allxdyn = zeros(nvar,length(tout),nsols(idp));
    allslope = zeros(nvar,length(tout),nsols(idp));
    
    % choose/determine points that have mss
    if isfield(allnss,sprintf('iid%d',idp));
        msspts = find(allnss.(['iid' num2str(idp)]));
        sslps = allnss.(['iid' num2str(idp)])(msspts);
        ss = unique(sslps);
        nss = length(ss);
        styles = {':','--','-.','-'};      
        heqfig = [];
        htfig = [];
        h3a = [];
        ha = [];
        
        for iss = 1:nss
            allmsspts = msspts(sslps==ss(iss));
            allisspvec = alliidpvec(allmsspts,:,idp);
            allissxeq = alliidxeq(:,allmsspts,idp);
            nmsspts = length(allmsspts);
            Line.LineStyle = styles{iss};
                                               
            % collect corresponding limit points
            for ipt = 1:nmsspts
                index =...
                cat(1,siid.(['iid' num2str(idp)]).(['pt' num2str(allmsspts(ipt))]).s1.index);
                s1 =...
                siid.(['iid' num2str(idp)]).(['pt' num2str(allmsspts(ipt))]).s1;
                x1ind =...
                siid.(['iid' num2str(idp)]).(['pt' num2str(allmsspts(ipt))]).x1(:,index);
                x1 =...
                siid.(['iid' num2str(idp)]).(['pt' num2str(allmsspts(ipt))]).x1;
                flux1 =...
                siid.(['iid' num2str(idp)]).(['pt' num2str(allmsspts(ipt))]).flux;
                eigind =...
                siid.(['iid' num2str(idp)]).(['pt' num2str(allmsspts(ipt))]).f1(:,index);
                f1 =...
                siid.(['iid' num2str(idp)]).(['pt' num2str(allmsspts(ipt))]).f1;
                xLPval = x1ind(1:nvar,:);
                pLPval = x1ind(nvar+1:end,:);
                pvec = allisspvec(1,:);
                
                % get saddle point
                eps1 = 1e-4;
                saddle = [];
                while isempty(saddle)
                    [saddle,saddlepar] = getsaddlenode(s1,x1,eps1);
                    eps1 = eps1*10;
                end
                
                % get eigen value and eigne vector at saddle
                pvec(ap) = saddlepar;
                model.PM(ac-length(saddle)) = saddlepar;
                [~,lambda,w] = getKotteJacobian(saddle,pvec,model);
                
                % get 2 steady states by perturbing initial states from
                % saddle
                SAxeq = [];
                iter = 1;
                while iter <= 2
                    delx = [1 -1;1 -1;1 -1];
                    ival = saddle+1e-2*delx(:,iter);
                    [~,xeq] = solveODEonly(1,ival,model,pvec,opts,tout);
                    if ~isempty(SAxeq)
                        if ~any(abs(SAxeq-repmat(xeq,1,size(SAxeq,2)))<=1e-8)
                            SAxeq = [SAxeq xeq];
                        end
                    else
                        SAxeq = xeq;
                    end   
                    if abs(xeq-SAxeq(:,1))<1e-8
                        Line.Color = colorSpec{1};
                        Point.MarkerEdgeColor = [0 0 0];
                        Point.MarkerFaceColor = colorSpec{1};
                    elseif size(SAxeq,2)>1 & abs(xeq-SAxeq(:,2))<1e-8
                        Line.Color = colorSpec{2};
                        Point.MarkerEdgeColor = [0 0 0];
                        Point.MarkerFaceColor = colorSpec{2};
                    else
                        Point = [];
                    end 
                    [heqfig,h3a] =...
                    FIGmssEqIvalPerturbations(ival,xeq,2,[1 2],heqfig,h3a,Point); 
                    iter = iter+1;
                end
                
                
                
                % separatrix curves through numerical approximation from
                % saddle
%                 NumericalSeparatrix(model,pvec,opts,ap,s1,x1,SAxeq);
                  
                
%                 test code for quiver3 - run on server - memory intensive
%                 abandon
%                 pvec(ap) = 0.1450;
%                 [X, Y, Z] = meshgrid(0:0.1:6, 0:0.1:2,0:0.1:4);
%                 givenModel = @(t,x)KotteODE(t,x,model,pvec);
%                 DX = zeros(size(X,1),size(X,2),size(X,3));
%                 DY = zeros(size(Y,1),size(Y,2),size(Y,3));
%                 DZ = zeros(size(Z,1),size(Z,2),size(Z,3));
%                 for i = 1:size(X,1)
%                     for j = 1:size(X,2)
%                         for k = 1:size(X,3)
%                             dX = givenModel(0,[X(i,j,k);Y(i,j,k);Z(i,j,k)]);
%                             DX(i,j,k) = dX(1);
%                             DY(i,j,k) = dX(2);
%                             DZ(i,j,k) = dX(3);
%                         end
%                     end
%                 end
%                 [gradX,gradY,gradZ] = gradient(DX,0.1,0.1,0.1);
                                
%                 bifurcationPlot(x1,x1(nvar+1:end,:),s1,f1,2,1);
                % choose a var index in xLPval that has a zero
                % /positive real eigen value
%                 [poseigind,~] = find(eigind>=0);
%                 poseigind = unique(poseigind);
                nindexpts = length(index);
                
                % choose nxsspts random points for poseeigind within bounds in
                % xLPval                
%                 for kindexpts = 1:(nindexpts-1)
%                     pxLPval =...
%                     sampleEKP(xLPval(:,kindexpts)',xLPval(:,kindexpts),...
%                               xLPval(:,kindexpts+1),[1 2 3],nxsspts)';
%                     ppLPval =...
%                     sampleEKP(pLPval(:,kindexpts)',pLPval(1,kindexpts),...
%                               pLPval(1,kindexpts+1),1,nxsspts);
                    
                    % subject these points to perturbations  
%                     for ixsspt = 1:nxsspts
%                         ival = pxLPval(:,ixsspt); 
%                         pvec(ap) = ppLPval(ixsspt);
%                         [xdyn,xeq,fdyn,feq,slope] =...
%                         solveODEonly(1,ival,model,pvec,opts,tout); 
%                         allival(:,ipt) = ival;
%                         allxeq(:,isol) = xeq;
%                         allxdyn(:,:,isol) = xdyn;
%                         allslope(:,:,isol) = slope;
%                         isol = isol+1;
%                         
%                         Line.DisplayName =...
%                         ['pt' num2str(ipt) 'smp' num2str(ixsspt)];       
%                     
%                         if abs(xeq-LPxeq(:,1))<1e-8
%                             Line.Color = colorSpec{1};
%                             Point.MarkerEdgeColor = [0 0 0];
%                             Point.MarkerFaceColor = colorSpec{1};
%                         elseif abs(xeq-LPxeq(:,2))<1e-8
%                             Line.Color = colorSpec{2};
%                             Point.MarkerEdgeColor = [0 0 0];
%                             Point.MarkerFaceColor = colorSpec{2};
%                         end  
                        
%                         [heqfig,h3a] =...
%                         FIGmssEqIvalPerturbations(ival,xeq,2,[],heqfig,h3a,Point); 
                        
                        % plot dynbamic trajectories
%                         [htfig,ha] =...
%                         FIGodetrajectories(xdyn,ival,xeq,2,[1 2],htfig,ha,[],Point);
                        
                        % Plot all dynamic information
%                         [hfig,hsfig,hline] =...
%                         FIGmssdynamicswpvec(xdyn,tout,pvec,[1 2 3],1,1,2,1:find(tout==100),Line,hfig,hsfig,hline);
%                     end
%                 end                
            end            
        end
        % print the equlibrium points in different color
        % plot3 returns line object handles        

%         for ixeq = 1:size(LPxeq,2)
%             plot3(h3a,LPxeq(1,ixeq),LPxeq(2,ixeq),LPxeq(3,ixeq),...
%                                               'Marker','o',...
%                                               'MarkerSize',10,...
%                                               'MarkerFaceColor','r',...
%                                               'MarkerEdgeColor','r');
%         end
% Figures from simulation due to changes in parameters - equilibrium and dynamic
% Plot trajectories and trialing sepratrix vector plots
    end
end
