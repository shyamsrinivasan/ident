% load('C:\Users\shyam\Documents\Courses\CHE1125Project\Results\KotteModel\RegulationVariation_July30.mat');
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
    for ialpha1 = 1:nalpha2
        startid = find(allpvec(:,idp(1))==alpha2_range(ialpha1),1,'first');
        endid = find(allpvec(:,idp(1))==alpha2_range(ialpha1),1,'last');
        if isfield(allnss,sprintf('iid%d',iid))
            % find all alpha1 that are bistable
            bistable_alpha1(1,ialpha1) =...
            max(allpvec(allnss.(['iid' num2str(iid)])(startid:endid)==4,idp(2)));
            bistable_alpha1(2,ialpha1) =...
            min(allpvec(allnss.(['iid' num2str(iid)])(startid:endid)==4,idp(2)));
            % find all alpha1 that are monostable with unstable regions
            unstable_alpha1(1,ialpha1) =...
            max(allpvec(allnss.(['iid' num2str(iid)])(startid:endid)==3,idp(2)));
            unstable_alpha1(2,ialpha1) =...
            min(allpvec(allnss.(['iid' num2str(iid)])(startid:endid)==3,idp(2)));
            % find all alpha1 that are monostable
            mstable_alpha1(1,ialpha1) =...
            max(allpvec(allnss.(['iid' num2str(iid)])(startid:endid)==0,idp(2)));
            mstable_alpha1(2,ialpha1) =...
            min(allpvec(allnss.(['iid' num2str(iid)])(startid:endid)==0,idp(2)));
        end
    end
end
figure
plot(log(repmat(alpha2_range',2,1)),log(mstable_alpha1),'LineStyle','none','Marker','.','MarkerSize',10);
hold on
plot(log(repmat(alpha2_range',2,1)),log(bistable_alpha1),'LineStyle','none','Marker','.','MarkerSize',10);
plot(log(repmat(alpha2_range',2,1)),log(unstable_alpha1),'LineStyle','none','Marker','.','MarkerSize',10);
xlabel('log(Lfbp)');
ylabel('log(KFbpPEP)');

figure
plot(repmat(alpha2_range',2,1),mstable_alpha1,'LineStyle','none','Marker','.','MarkerSize',10);
hold on
plot(repmat(alpha2_range',2,1),bistable_alpha1,'LineStyle','none','Marker','.','MarkerSize',10);
plot(repmat(alpha2_range',2,1),unstable_alpha1,'LineStyle','none','Marker','.','MarkerSize',10);
xlabel('log(Lfbp)');
ylabel('log(KFbpPEP)');

% bifurcation diagrams for mss values only
for iid = 1:ndp
    allpvec = alliidpvec(:,:,iid);
    % unique alpha2
    alpha2_range = unique(allpvec(:,idp(1)));
    nalpha2 = length(alpha2_range);
    % max and min bistable apha1 for every alpha2
    bistable_alpha1 = zeros(2,nalpha2);
    unstable_alpha1 = zeros(2,nalpha2);
    mstable_alpha1 = zeros(2,nalpha2);
    last_endid = 0;
    for ialpha1 = 1:nalpha2
        hf1 = [];
        startid = find(allpvec(:,idp(1))==alpha2_range(ialpha1),1,'first');
        endid = find(allpvec(:,idp(1))==alpha2_range(ialpha1),1,'last');
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
        last_endid = endid;
    end
end