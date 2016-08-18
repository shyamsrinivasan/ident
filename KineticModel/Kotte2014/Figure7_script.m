% data from Figure 5
load('C:\Users\shyam\Documents\Courses\CHE1125Project\Results\KotteModel\kPEPoutVariation_Aug17.mat');

% needed variables: alliidpvec,alliidxeq,alliidfeq,tout,ap,idp;
npts = size(alliidpvec,1);
nvar = size(alliidxeq,1);
ndp = size(alliidpvec,3);
tout = tspan;
% umber of bifurcation points
nbifpts = 6000;

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
            for ipt = 1:npts
                addanot.text = ['R' num2str(ipt)];
                % if point capable of mss
                if ismember(ipt,allmsspts)   
%                     s1 =...
%                     siid.(['iid' num2str(iid)]).(['pt' num2str(ipt)]).s1;
%                     x1 =...
%                     siid.(['iid' num2str(iid)]).(['pt' num2str(ipt)]).x1;
%                     f1 =...
%                     siid.(['iid' num2str(iid)]).(['pt' num2str(ipt)]).f1;
%                     index =...
%                     cat(1,siid.(['iid' num2str(iid)]).(['pt' num2str(ipt)]).s1.index);
                    allacetate(ipt,:) =...
                    siid.(['iid' num2str(iid)]).(['pt' num2str(ipt)]).x1(end,:);
%                     hf1 = bifurcationPlot(x1,s1,f1,[4,1],ap,hf1,addanot);
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
        end
    end
end