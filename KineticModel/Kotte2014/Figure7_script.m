% data from Figure 5
load('C:\Users\shyam\Documents\Courses\CHE1125Project\Results\KotteModel\kPEPoutVariation_Aug19.mat');

% needed variables: alliidpvec,alliidxeq,alliidfeq,tout,ap,idp,(allacxeq,allacfeq,acetate);
npts = size(alliidpvec,1);
nvar = size(alliidxeq,1);
nflux = size(alliidfeq,1);
ndp = size(alliidpvec,3);
tout = tspan;
tspanf = 0:0.1:2000;

% umber of continuation points
nbifpts = 6000;
ac = find(strcmpi(model.mets,'ac[e]'));

for iid = 1:ndp
    if isfield(allnss,sprintf('iid%d',iid))
        msspts = find(allnss.(['iid' num2str(iid)]));
        sslps = allnss.(['iid' num2str(iid)])(msspts);
        ss = unique(sslps);
        nss = length(ss);
        allmsspts = [];
        for iss = 1:nss
            hf1 = [];
            allmsspts = union(allmsspts,msspts(sslps==ss(iss)));
            allacetate = zeros(length(allmsspts),nbifpts);
            bistablex = zeros(2*nvar,npts,length(acetate));
            bistablef = zeros(2*nflux,npts,length(acetate));
            % unistablex = zeros(nvar,npts,length(acetate));
            % unistablef = zeros(nflux,npts,length(acetate));
            saddleac = zeros(npts,length(acetate));
            for ipt = 1:npts
                pvec = alliidpvec(ipt,:,iid);
                addanot.text = ['R' num2str(ipt)];
                % if point capable of mss
                if ismember(ipt,allmsspts)  
                    % find two steady states for different acetate
                    % concentrations
%                     s1 =...
%                     siid.(['iid' num2str(iid)]).(['pt' num2str(ipt)]).s1;
%                     x1 =...
%                     siid.(['iid' num2str(iid)]).(['pt' num2str(ipt)]).x1;
                    s1.pt1 =...
                    siid.(['iid' num2str(iid)]).(['pt' num2str(ipt)]);
                    for iac = 1:length(acetate)
                        fprintf('%d\n',iac); % debug
                        % calculate saddle for each acetate concentration
                        eps = 1e-4;
                        saddle = [];
                        saddlepar = [];
                        while isempty(saddlepar)
                            [saddle,saddlepar,status] = geteqcontpt(s1,acetate(iac),eps);
                            eps = eps*10;
                        end
                        if status % obtained good saddle node
                            pvec(ap) = saddlepar;
                            saddleac(ipt,iac) = saddlepar;
                            model.PM(ac-length(saddle)) = saddlepar;
                            % perturb saddle to get steady states
                            eps = 1e-4;                            
                            pival = saddle+eps*[1;1;1];
                            [~,xeq1,~,feq1] =...
                            solveODEonly(1,pival,model,pvec,opts,tspanf);
                            nival = saddle-eps*[1;1;1];
                            [~,xeq2,~,feq2] =...
                            solveODEonly(1,nival,model,pvec,opts,tspanf);
                            bistablex(1:nvar,ipt,iac) = xeq1;
                            bistablex(nvar+1:end,ipt,iac) = xeq2;
                            bistablef(1:nflux,ipt,iac) = feq1;
                            bistablef(nflux+1:end,ipt,iac) = feq2;
                        elseif ~status 
                            % saddle parameter out of bifurcation bounds
                            % get the one possible steady state
                            pvec(ap) = acetate(iac);
                            saddleac(ipt,iac) = acetate(iac);
                            model.PM(ac-length(saddle)) = acetate(iac);
                            [~,xeq1,~,feq1] =...
                            solveODEonly(1,M,model,pvec,opts,tspanf);
                            bistablex(1:nvar,ipt,iac) = xeq1;
                            bistablex(nvar+1:end,ipt,iac) = xeq1;
                            bistablef(1:nflux,ipt,iac) = feq1;
                            bistablef(nflux+1:end,ipt,iac) = feq1;
                        end
                    end
                    allacetate(ipt,:) =...
                    siid.(['iid' num2str(iid)]).(['pt' num2str(ipt)]).x1(end,:);
                % else if point not capable of mss
                elseif ~ismember(ipt,allmsspts)
                    % find the monostable steady states for different
                    % acetate concentrations
                    for iac = 1:length(acetate)
                        % take data from existing information in workspace
                        saddleac(ipt,iac) = acetate(iac);
                        bistablex(1:nvar,ipt,iac) = allacxeq(:,ipt,iac);
                        bistablex(nvar+1:end,ipt,iac) = allacxeq(:,ipt,iac);
                        bistablef(1:nflux,ipt,iac) = allacfeq(:,ipt,iac);
                        bistablef(nflux+1:end,ipt,iac) = allacfeq(:,ipt,iac);
                    end
                end

            end  
            endacetate = min(max(allacetate,[],2));
            % find indices for every row in allacetate close to endacetate
            [rwid,ind] = find(abs(allacetate-repmat(endacetate,length(allmsspts),nbifpts))<=1e-2);
            ind = [rwid ind];
            ind = sortrows(ind,1);
            [~,ia] = unique(ind(:,1),'rows','first');
            C = ind(ia,:);
            minac = zeros(length(allmsspts),1);
            maxac = zeros(length(allmsspts),1);   
            kpepval = zeros(length(allmsspts),1);
            for ipt = 1:npts
                if ismember(ipt,allmsspts)
                    x1 =...
                    siid.(['iid' num2str(iid)]).(['pt' num2str(ipt)]).x1(:,1:C(C(:,1)==ipt,2));
                    flux1 =...
                    siid.(['iid' num2str(iid)]).(['pt' num2str(ipt)]).flux(:,1:C(C(:,1)==ipt,2));
                    index =...
                    cat(1,siid.(['iid' num2str(iid)]).(['pt' num2str(ipt)]).s1.index);
                    ac1 = x1(end,:);
                    lpsval = x1(:,index(2):index(end-1));
                    minac(ipt==allmsspts) = min(lpsval(end,:));
                    maxac(ipt==allmsspts) = max(lpsval(end,:));
                    kpepval(ipt==allmsspts) = alliidpvec(ipt,idp,iid);
    %                 selectpts = x1(
                end
            end
            figure
            plot(minac,kpepval,'Color','r','LineWidth',2);
            hold on
            plot(maxac,kpepval,'Color','r','LineWidth',2);
            plot([max(minac) max(maxac)],[max(kpepval) max(kpepval)],...
                 'Color','r','LineWidth',2);
            xlabel('acetate a.u');
            ylabel('kPEPout a.u');
        end
    end
end

%% Plot kPEPout vs acetate for region of bistability
% load data
load('C:\Users\shyam\Documents\Courses\CHE1125Project\Results\KotteModel\kPEPoutVariation_Aug20.mat');

npts = size(alliidpvec,1);
nvar = size(alliidxeq,1);
nflux = size(alliidfeq,1);
ndp = size(alliidpvec,3);
tout = tspan;
tspanf = 0:0.1:2000;

figure
plot(minac,kpepval,'Color','r','LineWidth',2);
hold on
plot(maxac,kpepval,'Color','r','LineWidth',2);
plot([max(minac) max(maxac)],[max(kpepval) max(kpepval)],...
     'Color','r','LineWidth',2);
xlabel('acetate a.u');
ylabel('kPEPout a.u');

%% Plots for kPEPout and acetate vs PEP and v4 flux
% load data from simulation of previous cell
load('C:\Users\shyam\Documents\Courses\CHE1125Project\Results\KotteModel\kPEPoutVariation_Aug20.mat');

colorSpec = chooseColors(5,{'Green','Purple','Red','Navy','HotPink'});

nact = length(acetate);
for iac = 1:length(acetate)
    hc = figure;
    plot(alliidpvec(:,11),reshape(bistablex(1,:,iac),1,npts),...
        'Color',colorSpec{1},'Marker','.','MarkerSize',15);
    hold on
    plot(alliidpvec(:,11),reshape(bistablex(4,:,iac),1,npts),...
        'Color',colorSpec{2},'Marker','.','MarkerSize',15);
    xlabel('kPEPout');
    ylabel('PEP a.u');
    path = 'C:\Users\shyam\Documents\Courses\CHE1125Project\Results\KotteModel';
    filename = [path '\PEPvkPEPout_acetate' num2str(iac)];
    saveas(hc,filename,'epsc');
%     savefig(hc,filename);
    close(hc);
end



for iac = 1:length(acetate)
    hf = figure;
    plot(alliidpvec(:,11),reshape(bistablef(5,:,iac),1,npts),...
        'Color',colorSpec{1},'Marker','.','MarkerSize',15);
    hold on
    plot(alliidpvec(:,11),reshape(bistablef(10,:,iac),1,npts),...
        'Color',colorSpec{2},'Marker','.','MarkerSize',15);
    xlabel('kPEPout');
    ylabel('v4');
    path = 'C:\Users\shyam\Documents\Courses\CHE1125Project\Results\KotteModel';
    filename = [path '\v4vkPEPout_acetate' num2str(iac)];
    saveas(hf,filename,'epsc');
%     savefig(hf,filename);
    close(hf);
end
