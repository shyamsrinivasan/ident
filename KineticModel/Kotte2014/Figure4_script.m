% load('C:\Users\shyam\Documents\Courses\CHE1125Project\Results\KotteModel\RegulationVariation_July26.mat');
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

