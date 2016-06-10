% runDynamicsims
% simulate from one or more limit points or within their boundaries 
% as initial conditions to test stability/dynamic behaviour of points
% load the relevant simulation dataset
load('C:\Users\shyam\Documents\Courses\CHE1125Project\Results\KotteModel\VmaxVariationAll_May28.mat');
Figures from simulation due to changes in parameters - equilibrium and dynamic

% needed variables: alliidpvec,alliidxeq,alliidfeq,tout,ap;
npts = size(alliidpvec,1);
nvar = size(alliidxeq,1);
ndp = size(alliidpvec,3);

FIGmssdynamicswpvec(alliidxdyn,tout,alliidpvec,1,1,1:1000,'conc',...
                    1:find(tout==100))

opts = odeset('RelTol',1e-12,'AbsTol',1e-10);

% determine the #parameters/combinations that have been changed
for idp = 1:ndp
    % determine actual # parameters that have been changed
    npar = length(alliidpvec(1,:,idp));
    if npts>1
        diffpar = find(alliidpvec(1,:,idp)~=alliidpvec(2,:,idp));
    end    
    hfig = [];
    hsfig = [];
    hline = [];
    colorSpec = chooseColors(4);
    allival = [];
    everyxeq = [];
    
    % choose/determine points that have mss
    if isfield(allnss,sprintf('iid%d',idp));
        msspts = find(allnss.(['iid' num2str(idp)]));
        sslps = allnss.(['iid' num2str(idp)])(msspts);
        ss = unique(sslps);
        nss = length(ss);
        styles = {':','--','-.','-'};      
        heqfig = figure;
        h3a = [];
        
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
                s1 = siid.(['iid' num2str(idp)]).(['pt' num2str(allmsspts(ipt))]).s1;
                x1ind =...
                siid.(['iid' num2str(idp)]).(['pt' num2str(allmsspts(ipt))]).x1(:,index);
                x1 = siid.(['iid' num2str(idp)]).(['pt' num2str(allmsspts(ipt))]).x1;
                flux1 = siid.(['iid' num2str(idp)]).(['pt' num2str(allmsspts(ipt))]).flux;
                eigind =...
                siid.(['iid' num2str(idp)]).(['pt' num2str(allmsspts(ipt))]).f1(:,index);
                f1 = siid.(['iid' num2str(idp)]).(['pt' num2str(allmsspts(ipt))]).f1;
                xLPval = x1ind(1:nvar,:);
                pLPval = x1ind(nvar+1:end,:);
                pvec = allisspvec(1,:);
                
                % find the 2 or more steady states from the LPs
                LPxeq = [];
                for it = 1:size(xLPval,2)
                    ival = xLPval(:,it); 
                    pvec(ap) = pLPval(it);
                    [~,xeq] =...
                    solveODEonly(1,ival,model,pvec,opts,tout); 
                    if ~isempty(LPxeq)
                        if ~any(abs(LPxeq-repmat(xeq,1,size(LPxeq,2)))<=1e-8)
                            LPxeq = [LPxeq xeq];
                        end
                    else
                        LPxeq = xeq;
                    end                
                end
                                
%                 bifurcationPlot(x1,x1(nvar+1:end,:),s1,f1,2,1);
                % choose a var index in xLPval that has a zero
                % /positive real eigen value
%                 [poseigind,~] = find(eigind>=0);
%                 poseigind = unique(poseigind);
                nindexpts = length(index);
                
                % choose nxsspts random points for poseeigind within bounds in
                % xLPval
                nxsspts = 2;
                for kindexpts = 1:(nindexpts-1)
                    pxLPval =...
                    sampleEKP(xLPval(:,kindexpts)',xLPval(:,kindexpts),...
                              xLPval(:,kindexpts+1),[1 2 3],nxsspts)';
                    ppLPval =...
                    sampleEKP(pLPval(:,kindexpts)',pLPval(1,kindexpts),...
                              pLPval(1,kindexpts+1),1,nxsspts);
                    
                    % subject these points to perturbations  
                    for ixsspt = 1:nxsspts
                        ival = pxLPval(:,ixsspt); 
                        pvec(ap) = ppLPval(ixsspt);
                        [allxdyn,allxeq,allfdyn,allfeq] =...
                        solveODEonly(1,ival,model,pvec,opts,tout); 
                        allival = [allival ival];
                        everyxeq = [everyxeq allxeq];
                        Line.DisplayName =...
                        ['pt' num2str(ipt) 'smp' num2str(ixsspt)];
                        [heqfig,h3a] =...
                        FIGmssEqIvalPerturbations(ival,allxeq,2,[],heqfig,h3a);                        
%                         h3fig = plot3([ival(1) allxeq(1)],...
%                                     [ival(2) allxeq(2)],...
%                                     [ival(3) allxeq(3)],...
%                                       'Marker','o',...
%                                       'MarkerSize',6,...
%                                       'MarkerFaceColor',[.49 1 .63],...
%                                       'MarkerEdgeColor','k');
%                         hold on                
                    
                        if abs(allxeq-LPxeq(:,1))<1e-8
                            Line.Color = colorSpec{1};
                        elseif abs(allxeq-LPxeq(:,2))<1e-8
                            Line.Color = colorSpec{2};
                        end  
                        
                        % Plot all dynamic information
                        [hfig,hsfig,hline] =...
                        FIGmssdynamicswpvec(allxdyn,tout,pvec,[1 2 3],1,1,2,1:find(tout==100),Line,hfig,hsfig,hline);
                    end
                end                
            end            
        end
        % print the equlibrium points in different color
        % plot3 returns line object handles        
        for ixeq = 1:size(LPxeq,2)
            plot3(h3a,LPxeq(1,ixeq),LPxeq(2,ixeq),LPxeq(3,ixeq),...
                                              'Marker','o',...
                                              'MarkerSize',10,...
                                              'MarkerFaceColor','r',...
                                              'MarkerEdgeColor','r');
        end
% Figures from simulation due to changes in parameters - equilibrium and dynamic
    end
end
