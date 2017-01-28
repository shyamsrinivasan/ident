load('C:\Users\shyam\Documents\Courses\CHE1125Project\Results\KotteModel\RegulationVariation_July30.mat');
% idp = [3 4] % sampled together in a mesh grid
npts = size(alliidpvec,1);
nvar = size(alliidxeq,1);
ndp = size(alliidpvec,3);

if exist('tspan','var')
    tout = tspan;
end

for iid = 1:ndp
    allpvec = alliidpvec(:,:,iid);
    % unique alpha2
    alpha2_range = unique(allpvec(:,idp(1)));
    nalpha2 = length(alpha2_range);
    % max and min bistable apha1 for every alpha2
    bistable_alpha1 = zeros(2,nalpha2);
    unstable_alpha1 = zeros(2,nalpha2);
    mstable_alpha1 = zeros(2,nalpha2);
    for ialpha2 = 1:nalpha2
        startid = find(allpvec(:,idp(1))==alpha2_range(ialpha2),1,'first');
        endid = find(allpvec(:,idp(1))==alpha2_range(ialpha2),1,'last');
        if isfield(allnss,sprintf('iid%d',iid))
            % find all alpha1 that are bistable
            bistable_alpha1(1,ialpha2) =...
            max(allpvec(allnss.(['iid' num2str(iid)])(startid:endid)==4,idp(2)));
            bistable_alpha1(2,ialpha2) =...
            min(allpvec(allnss.(['iid' num2str(iid)])(startid:endid)==4,idp(2)));
            % find all alpha1 that are monostable with unstable regions
            unstable_alpha1(1,ialpha2) =...
            max(allpvec(allnss.(['iid' num2str(iid)])(startid:endid)==3,idp(2)));
            unstable_alpha1(2,ialpha2) =...
            min(allpvec(allnss.(['iid' num2str(iid)])(startid:endid)==3,idp(2)));
            % find all alpha1 that are monostable
            mstable_alpha1(1,ialpha2) =...
            max(allpvec(allnss.(['iid' num2str(iid)])(startid:endid)==0,idp(2)));
            mstable_alpha1(2,ialpha2) =...
            min(allpvec(allnss.(['iid' num2str(iid)])(startid:endid)==0,idp(2)));
        end
    end
end
% figure
% plot(log(repmat(alpha2_range',2,1)),log(mstable_alpha1),'LineStyle','none','Marker','.','MarkerSize',10);
% hold on
% plot(log(repmat(alpha2_range',2,1)),log(bistable_alpha1),'LineStyle','none','Marker','.','MarkerSize',10);
% plot(log(repmat(alpha2_range',2,1)),log(unstable_alpha1),'LineStyle','none','Marker','.','MarkerSize',10);
% xlabel('log(Lfbp)');
% ylabel('log(KFbpPEP)');

figure
plot(log(alpha2_range'),mstable_alpha1(1,:),'Color','k','LineWidth',2);
hold on
plot(log(alpha2_range'),mstable_alpha1(2,:),'Color','k','LineWidth',2);
plot(log(alpha2_range'),bistable_alpha1(1,:),'Color','r','LineWidth',2);
plot(log(alpha2_range'),bistable_alpha1(2,:),'Color','r','LineWidth',2);
plot(log(alpha2_range'),unstable_alpha1(1,:),'Color','b','LineWidth',2);
plot(log(alpha2_range'),unstable_alpha1(2,:),'Color','b','LineWidth',2);
xlabel('log(Lfbp)');
ylabel('KFbpPEP');

% pick one Lfbpo point and plot bifurcation diagrams for all corresponding
% KFbpPEP values
jalpha2 = size(alpha2_range,1);
if rem(jalpha2,2)~=0
    jalpha2 = (jalpha2+1)/2;
end

% bifurcation on KfbpPEP 
ap = 4;
fluxg = Kotte_givenFlux([M;model.PM],pvec,model);
for iid = 1:ndp
    allpvec = alliidpvec(:,:,iid);
    allxeq = alliidxeq(:,:,iid);
    % unique alpha2
    alpha2_range = unique(allpvec(:,idp(1)));
    nalpha2 = length(alpha2_range);    
    for ialpha2 = 1:nalpha2
        hf1 = [];
        pvec = allpvec(ialpha2,:);
        startid = find(allpvec(:,idp(1))==alpha2_range(ialpha2),1,'first');
        endid = find(allpvec(:,idp(1))==alpha2_range(ialpha2),1,'last');
        if ialpha2 == jalpha2
            Lfbp = alpha2_range(ialpha2);
            pvec(idp(1)) = alpha2_range(ialpha2);
            alpha1_range = unique(allpvec(startid:endid,idp(2)));
            nalpha1 = length(alpha1_range);   
            ialpha1 = 0;
            while (startid+ialpha1) < endid
                pvec(idp(2)) = alpha1_range(ialpha1+1);
                xeq = allxeq(:,startid+ialpha1);
                addanot.text = ['R' num2str(ialpha1+1)];
                % continue of new parameter
                ap = idp(2);
                 % run MATCONT
                [data,y,p] = execMATCONT(xeq,pvec,ap,fluxg,model);
                if ~isempty(data) && size(data.s1,1)>2
                    hf1 = bifurcationPlot(data.x1,data.s1,data.f1,[4,1],ap,hf1,addanot);
                end
                ialpha1 = ialpha1+1;
            end       
        end              
    end
end

% bifurcation on acetate diagrams for mss values only
for iid = 1:ndp
    allpvec = alliidpvec(:,:,iid);
    % unique alpha2
    alpha2_range = unique(allpvec(:,idp(1)));
    nalpha2 = length(alpha2_range);
    % max and min bistable apha1 for every alpha2
%     bistable_alpha1 = zeros(2,nalpha2);
%     unstable_alpha1 = zeros(2,nalpha2);
%     mstable_alpha1 = zeros(2,nalpha2);
    last_endid = 0;
    for ialpha2 = 1:nalpha2
        hf1 = [];
        startid = find(allpvec(:,idp(1))==alpha2_range(ialpha2),1,'first');
        endid = find(allpvec(:,idp(1))==alpha2_range(ialpha2),1,'last');
        if ialpha2 == jalpha2
            Lfbp = alpha2_range(ialpha2);
            if isfield(allnss,sprintf('iid%d',iid))
                mssid = find(allnss.(['iid' num2str(iid)])(startid:endid)==4);
                mssid = mssid+last_endid;
                for imss = 1:length(mssid)
                    data = siid.(['iid' num2str(iid)]).(['pt' num2str(mssid(imss))]);
    %                 s1 = siid.(['iid' num2str(iid)]).(['pt' num2str(mssid(imss))]).s1;
    %                 x1 = siid.(['iid' num2str(iid)]).(['pt' num2str(mssid(imss))]).x1;
    %                 f1 = siid.(['iid' num2str(iid)]).(['pt' num2str(mssid(imss))]).f1;
                    hf1 = bifurcationPlot(data.x1,data.s1,data.f1,[4,1],ap,hf1);
                end
            end
        end
        last_endid = endid;        
    end
end

%% get pep and v4 vs Kallos for different acetate 
acetate = orig_saddlepar;
ap = 9;
colorSpec = chooseColors(5,{'Green','Purple','Red','Navy','HotPink'});
saddleac = zeros(npts,length(acetate));
xeqac = zeros(2*nvar,npts,length(acetate));
feqac = zeros(2*length(fluxg),npts,length(acetate));
for iac = 1:length(acetate)    
    % calculate saddle for each acetate concentration
    eps = 1e-4;
    [saddle,saddlepar,status] = eqptwrapper(s,nvar,acetate(iac),eps);
    
    % good saddle node points only
    goodsaddle = saddle(:,logical(status));
    goodsaddlepar = saddlepar(logical(status));
    saddleac(logical(status),iac) = goodsaddlepar;
    allpvec(logical(status),ap) = goodsaddlepar;
    % saddle parameter out of bifurcation bounds
    % get the one possible steady state
    oubsaddle = saddle(:,~logical(status));
    oubsaddlepar = saddlepar(~logical(status));
    saddleac(~logical(status),iac) = acetate(iac);
    allpvec(~logical(status),ap) = acetate(iac);
    
    for ipt = 1:npts
        if ismember(ipt,find(status))
            model.PM(ac-length(orig_saddle)) = saddlepar(ipt); 
            % perturb saddle to get steady states
            eps = 1e-4;                            
            pival = saddle(:,ipt)+eps*[1;1;1];
            [~,xeq1,~,feq1] =...
            solveODEonly(1,pival,model,allpvec(ipt,:),opts,tspanf);
            nival = saddle(:,ipt)-eps*[1;1;1];
            [~,xeq2,~,feq2] =...
            solveODEonly(1,nival,model,allpvec(ipt,:),opts,tspanf);            
            xeqac(nvar+1:end,ipt,iac) = xeq2;
            feqac(length(fluxg)+1:end,ipt,iac) = feq2;
        else
            % get the only possible steady state
            model.PM(ac-length(orig_saddle)) = acetate(iac);
            [~,xeq1,~,feq1] = solveODEonly(1,M,model,allpvec(ipt,:),opts,tspan);            
            xeqac(nvar+1:end,ipt,iac) = xeq1;
            feqac(length(fluxg)+1:end,ipt,iac) = feq1;
        end
        xeqac(1:nvar,ipt,iac) = xeq1;
        feqac(1:length(fluxg),ipt,iac) = feq1;
    end
    
    hc1 = figure;
    hold on
    plot(allpvec(:,idp),xeqac(1,:),'Color',colorSpec{1},'LineWidth',2);
    plot(allpvec(:,idp),xeqac(4,:),'Color',colorSpec{2},'LineWidth',2);
    xlabel('k4cat s-1');
    ylabel('PEP a.u.');
    
    hc2 = figure;
    hold on
    plot(allpvec(:,idp),feqac(5,:),'Color',colorSpec{1},'LineWidth',2);
    plot(allpvec(:,idp),feqac(10,:),'Color',colorSpec{2},'LineWidth',2);
    xlabel('k4cat s-1');
    ylabel('v4 a.u.');
end