% plot figures from sims data
% load('C:\Users\shyam\Documents\Courses\CHE1125Project\Results\KotteModel\VmaxVariation_Jun01.mat');
% load('C:\Users\shyam\Documents\Courses\CHE1125Project\Results\KotteModel\VmaxVariation_Jun01_PerturbationSims_Jun17.mat');
npts = size(alliidpvec,1);
nvar = size(alliidxeq,1);
ndp = 1;
colorSpec = chooseColors(4,{'Green','Purple','Red','Orange'});
opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
ac = find(strcmpi(model.mets,'ac[e]'));

% plot trajectory figures
nsols = size(allival,2);

for idp = 1:ndp
    % determine actual # parameters that have been changed
    npar = length(alliidpvec(1,:,idp));
    if npts>1
        diffpar = find(alliidpvec(1,:,idp)~=alliidpvec(2,:,idp));
    end    
    hfig = [];
    hsfig = [];
    hline = [];    
    isol = 1;
    
    % choose/determine points that have mss
    if isfield(allnss,sprintf('iid%d',idp));
        msspts = find(allnss.(['iid' num2str(idp)]));
        sslps = allnss.(['iid' num2str(idp)])(msspts);
        ss = unique(sslps);
        nss = length(ss);
        styles = {':','--','-.','-'};      
        heqfig = figure;
        heqfig2 = figure;
        heqfig3 = [];
        htfig = [];
        htfig2 = [];
        htfig3 = [];
        h3a = [];
        h3a2 = [];
        h3a3 = [];
        ha = [];
        ha2 = [];
        ha3 = [];
        
        % calculation of all possible equilibitum points
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
                
                % find the 2 or more steady states from the LPs
                LPxeq = [];
                for it = 1:size(xLPval,2)
                    ival = xLPval(:,it); 
                    pvec(ap) = pLPval(it);
%                     model.PM(ac-length(ival)) = pLPval(it);
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
                bifurcationPlot(x1,x1,s1,f1,2,1,7);
            end              
        end
        for isol = 1:150
            ival = allival(:,isol);
            xeq = allxeq(:,isol);
            xdyn = allxdyn(:,:,isol);
            slope = allslope(:,:,isol);

            if abs(xeq-LPxeq(:,1))<1e-8
                Line.Color = colorSpec{1};
                Point.MarkerEdgeColor = [0 0 0];
                Point.MarkerFaceColor = colorSpec{1};
            elseif abs(xeq-LPxeq(:,2))<1e-8
                Line.Color = colorSpec{2};
                Point.MarkerEdgeColor = [0 0 0];
                Point.MarkerFaceColor = colorSpec{2};
            end 

            % plot equilibrium points and initial values in 2D or 3D
            [heqfig,h3a] =...
            FIGmssEqIvalPerturbations(ival,xeq,2,[1 2],heqfig,h3a,Point); 

            [heqfig2,h3a2] =...
            FIGmssEqIvalPerturbations(ival,xeq,2,[2 3],heqfig2,h3a2,Point); 

            [heqfig3,h3a3] =...
            FIGmssEqIvalPerturbations(ival,xeq,2,[1 3],heqfig3,h3a3,Point); 

            % plot dynamic trajectories in 2D or 3D
            [htfig,ha] =...
            FIGodetrajectories(xdyn,ival,xeq,2,[1 2],htfig,ha,Line,Point);

            [htfig2,ha2] =...
            FIGodetrajectories(xdyn,ival,xeq,2,[2 3],htfig2,ha2,Line,Point);

            [htfig3,ha3] =...
            FIGodetrajectories(xdyn,ival,xeq,2,[1 3],htfig3,ha3,Line,Point);
        end
    end
end

