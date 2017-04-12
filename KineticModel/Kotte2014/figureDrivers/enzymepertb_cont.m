% bifurcation studies for enzyme perturbation cases using MATCONT
% this code tests whether the ability of a model to be bistable is governed
% by its bifucation w.r.t enzyme parameter
% values
runKotte
% get saddle node to get both stable steady states and get a bistable
% acetate concentration for perturbation
[orig_saddle,orig_saddlepar] = getsaddlenode(data.s1,data.x1,5e-3);
pvec(ap) = orig_saddlepar;
model.PM(ac-length(orig_saddle)) = orig_saddlepar;

% perturb saddle to get steady states
eps = 1e-4;
tspanf = 0:0.1:2000;
pival = orig_saddle+eps*[1;1;1];
[~,xeq1,~,feq1] = solveODEonly(1,pival,model,pvec,opts,tspanf);
nival = orig_saddle-eps*[1;1;1];
[~,xeq2,~,feq2] = solveODEonly(1,nival,model,pvec,opts,tspanf);
xss = [xeq1 xeq2];
%% continuation on enzyme parameters at original parameters (acetate & enzyme)
idp = [12 13 14];
pmin = zeros(length(idp),npts); % min
pmax = zeros(length(idp),npts); % max
origpvec1 = pvec;
figure
for id = 1:length(idp)    
    % xbounds = zeros(nvar,npts);
    ap = idp(id);
    
    % fix parameter at a lower value than 1
    if ap~=14
        pvec(ap) = 0.01;
    else
        pvec(ap) = 2;
    end
    
    % find equilibrium solution and continue
    [~,id1xeq] =...
    solveODEonly(npts,M,model,pvec,opts,tspan);
          
    % continute using MATCONT
    [s,mssid,nss] = setupMATCONT(@KotteMATCONT,@Kottecont_fluxcalc,...
                                @getKotteaxislabels,id1xeq,pvec,ap,model,...
                                fluxg,npts,1500,1);
    
    % identify boundaries of bistable enzyme parameters
    for ipt = 1:npts
        if ismember(ipt,mssid)
            index = cat(1,s.(['pt' num2str(ipt)]).s1.index);
            x1 = s.(['pt' num2str(ipt)]).x1;
            pcont = x1(nvar+1:end,index);
            pmin(id,ipt) = min(pcont(:,2:end-1));
            pmax(id,ipt) = max(pcont(:,2:end-1));
        end
    end    
    % reset parameters 
    pvec = origpvec1;
end
%% continuation on enzyme parameters at saddle node for nominal enzyme 
origpvec = pvec; % pvec backup
npts = 1;
pminsaddle = zeros(length(idp),npts); % min
pmaxsaddle = zeros(length(idp),npts); % max
for id = 1:length(idp)    
    ap = idp(id);
    
    % fix parameter at a lower value than 1
    if ap~=14
        pvec(ap) = 0.01;
    else
        pvec(ap) = 5;
    end
    
    % find equilibrium solution and continue
    [~,id2xeq] =...
    solveODEonly(npts,M,model,pvec,opts,tspan);
          
    % continute using MATCONT
    [s,mssid,nss] = setupMATCONT(@KotteMATCONT,@Kottecont_fluxcalc,...
                                @getKotteaxislabels,id2xeq,pvec,ap,model,...
                                fluxg,npts,1500,1);
    pvec = origpvec;
    
    % determine bistable parameter boundaries
    for ipt = 1:npts
        if ismember(ipt,mssid)
            index = cat(1,s.(['pt' num2str(ipt)]).s1.index);
            x1 = s.(['pt' num2str(ipt)]).x1;
            pcont = x1(nvar+1:end,index);
            pminsaddle(id,ipt) = min(pcont(:,2:end-1));
            pmaxsaddle(id,ipt) = max(pcont(:,2:end-1));
        end
    end   
end
% reset ap
ap = 9;
%% continuation on enzyme parameters at saddle node for enzyme perturbations
% sample parameters indicated by indices in idp
% cmb = [.05 1 1;1 .05 1;1 1 .05;.05 .05 .05;...
%        .125 1 1;1 .125 1;1 1 .125;.125 .125 .125;...
%        .25 1 1;1 .25 1;1 1 .25;.25 .25 .25;...
%        .5 1 1;1 .5 1;1 1 .5;.5 .5 .5;...
%        2 1 1;1 2 1;1 1 2;2 2 2;...
%        4 1 1;1 4 1;1 1 4;4 4 4];
% 
% type = 'together';
% npts = size(cmb,1);
% 
% if strcmpi(type,'together')
%     alliidpvec = zeros(npts,length(pvec),size(idp,1));
%     alliidxeq = zeros(length(M),npts,size(idp,1));
%     alliidxdyn = zeros(length(M),length(tspan),npts,size(idp,1));
%     alliidfeq = zeros(length(fluxg),npts,size(idp,1));
%     alliidfdyn = zeros(length(fluxg),length(tspan),npts,size(idp,1));
% else    
%     alliidpvec = zeros(npts,length(pvec),length(idp));
%     alliidxeq = zeros(length(M),npts,length(idp));
%     alliidxdyn = zeros(length(M),length(tspan),npts,length(idp));
%     alliidfeq = zeros(length(fluxg),npts,length(idp));
%     alliidfdyn = zeros(length(fluxg),length(tspan),npts,length(idp));
% end
% 
% % set acetate conentration
% pvec = [KEacetate,KFbpFBP,Lfbp,KFbpPEP,...
%         KEXPEP,vemax,KeFBP,ne,acetate,d,...
%         kPEPout,kEcat,vFbpmax,vEXmax];
% pvec(ap) = orig_saddlepar;
% model.PM(ac-length(orig_saddle)) = orig_saddlepar;
% 
% % set parameters from cmb at idp position(s)
% allpvec = repmat(pvec,npts,1);
% allpvec(:,idp) = cmb;
% allpvec(:,ap) = pvec(ap);
% 
% % continuation on enzyme parameters 12, 13 and 14 for perturbed values in
% % cmb
% % change ipt from 1 through 4 to cycle through different types of perturbations
% % see Table 1 in Methods for the types fo perturbations
% ipt = 4;
% hf1 = [];
% hf2 = [];
% hf3 = [];
% orig2pvec = allpvec;
% npts = 1;
% while ipt<size(cmb,1)
%     pvec = allpvec(ipt,:);
%     id2xeq = zeros(length(M),npts);    
%     id2feq = zeros(length(fluxg),npts);
%     addanot.text = ['E' num2str(ipt)];
%     for ip = 1:length(idp)
%         ap = idp(ip);   
%         if ap~=14
%             pvec(ap) = 0.01;
%         else
%             pvec(ap) = 5;
%         end
%         % get equilibrium point
%         [~,id2xeq,~,id2feq] =...
%         solveODEonly(npts,M,model,pvec,opts,tspan,...
%               [],id2xeq,[],id2feq);
%         % run MATCONT
%         [data,y,p] = execMATCONT(id2xeq,pvec,ap,fluxg,model);
%         if ~isempty(data) && size(data.s1,1)>2
%             if ap == 12
%                 hf1 = bifurcationPlot(data.x1,data.s1,data.f1,[4,1],ap,hf1,addanot);
%             elseif ap == 13
%                 hf2 = bifurcationPlot(data.x1,data.s1,data.f1,[4,1],ap,hf2,addanot);
%             elseif ap == 14
%                 hf3 = bifurcationPlot(data.x1,data.s1,data.f1,[4,1],ap,hf3,addanot);
%             end            
%         end
%         % save MATCONT results
%         s.(['pt' num2str(ipt)]) = data;
%         pvec = allpvec(ipt,:);
%     end    
%     ipt = ipt+4;
% end

%% continuation on enzyme parameters for perturbed values in cmb at saddle 
% node (same as above cell)
cmb = [.05 1 1;1 .05 1;1 1 .05;.05 .05 .05;...
       .125 1 1;1 .125 1;1 1 .125;.125 .125 .125;...
       .25 1 1;1 .25 1;1 1 .25;.25 .25 .25;...
       .5 1 1;1 .5 1;1 1 .5;.5 .5 .5;...
       2 1 1;1 2 1;1 1 2;2 2 2;...
       4 1 1;1 4 1;1 1 4;4 4 4];

type = 'together';
npts = size(cmb,1);

% set acetate conentration
k1cat = 1;
K1ac = 0.1;    % or 0.02
K3fdp = 0.1;
v3max = 1;
L3 = 4e6;
K3pep = 0.1;
v2max = 1;
K2pep = 0.3;
vemax = 1.1;        % for bifurcation analysis: 0.7:0.1:1.3
KeFDP = 0.45;       % or 0.45
ne = 2;             % or 2
acetate = 0.1;      % a.u acetate
d = 0.25;           % or 0.25 or 0.35
k4cat = 0.2;
pvec = [K1ac,K3fdp,L3,K3pep,...
        K2pep,vemax,KeFDP,ne,acetate,d,...
        k4cat,k1cat,v3max,v2max];
pvec(ap) = orig_saddlepar;
model.PM(ac-length(orig_saddle)) = orig_saddlepar;

% set parameters from cmb at idp position(s)
allpvec = repmat(pvec,npts,1);
allpvec(:,idp) = cmb;
allpvec(:,ap) = pvec(ap);

orig2pvec = allpvec;
for ip = 1:length(idp)
    hbif1 = [];
    ap = idp(ip);    
    
    % fix parameter at a lower value than 1
    if ap~=14
        allpvec(:,ap) = 0.01;
    else
        allpvec(:,ap) = 5;
    end
    
    % find equilibirum point for all npts points 
    [~,id3xeq] = solveODEonly(npts,M,model,allpvec,opts,tspan);

    % continue using MATCONT for all npts points
    [s,mssid,nss,hbif1] = setupMATCONT(@KotteMATCONT,@Kottecont_fluxcalc,...
                            @getKotteaxislabels,id3xeq,allpvec,ap,model,...
                            fluxg,npts,900,1,hbif1); 
    
    % store continuation data
    cmbsiid.(['iid' num2str(ip)]) = s;
    cmbmssid.(['iid' num2str(ip)]) = mssid;
    % allnss.(['iid' num2str(ip)]) = nss;
    
     % reset parameter
    allpvec = orig2pvec;
end
% reset ap
ap = 9;
%% identify bistable bounds
fprintf('Summary of Enzyme Parameter Continuation\n');
for ip = 1:length(fieldnames(cmbsiid))
    fprintf('\nContinuation on parameter %d\n:',idp(ip));
    fprintf('Bistable Perturbation Vectors :\n');
    mssid = cmbmssid.(['iid' num2str(ip)]);    
    pbounds = zeros(2,length(mssid));
    xbounds = zeros(nvar,2*length(mssid));
    mssipt = 1;
    for ipt = 1:npts
        if ismember(ipt,mssid)
            index = cat(1,cmbsiid.(['iid' num2str(ip)]).(['pt' num2str(ipt)]).s1.index);
            x1 = cmbsiid.(['iid' num2str(ip)]).(['pt' num2str(ipt)]).x1;
            pcont = x1(nvar+1:end,index);
            pbounds(1,mssipt) = min(pcont(:,2:end-1));
            pbounds(2,mssipt) = max(pcont(:,2:end-1));
            fprintf('%d\t',ipt);
            mssipt = mssipt + 1;
        end
    end
    pbistable.(['iid' num2str(ip)]) = pbounds;    
end

%% continuation on acetate for perturbed enzyme parameter values
cmb = [.05 1 1;1 .05 1;1 1 .05;.05 .05 .05;...
       .125 1 1;1 .125 1;1 1 .125;.125 .125 .125;...
       .25 1 1;1 .25 1;1 1 .25;.25 .25 .25;...
       .5 1 1;1 .5 1;1 1 .5;.5 .5 .5;...
       2 1 1;1 2 1;1 1 2;2 2 2;...
       4 1 1;1 4 1;1 1 4;4 4 4];
type = 'together';
npts = size(cmb,1);

% set acetate conentration at lowest value for continuation after
k1cat = 1;
K1ac = 0.1;    % or 0.02
K3fdp = 0.1;
v3max = 1;
L3 = 4e6;
K3pep = 0.1;
v2max = 1;
K2pep = 0.3;
vemax = 1.1;        % for bifurcation analysis: 0.7:0.1:1.3
KeFDP = 0.45;       % or 0.45
ne = 2;             % or 2
acetate = 0.1;      % a.u acetate
d = 0.25;           % or 0.25 or 0.35
k4cat = 0.2;
pvec = [K1ac,K3fdp,L3,K3pep,...
        K2pep,vemax,KeFDP,ne,acetate,d,...
        k4cat,k1cat,v3max,v2max];
pvec(ap) = 0.1;
model.PM(ac-length(orig_saddle)) = 0.1;

% set parameters from cmb at idp position(s)
allpvec = repmat(pvec,npts,1);
allpvec(:,idp) = cmb;
allpvec(:,ap) = pvec(ap);

orig2pvec = allpvec;
npts = 1;
ipt = 4;
while ipt<=size(cmb,1)
% for ip = 1:npts
    hbif2 = figure;
    
%         if ap~=14
%             pvec(ap) = 0.01;
%         else
%             pvec(ap) = 5;
%         end

    % find equilibirum point for all npts points 
    [~,id4xeq] = solveODEonly(npts,M,model,allpvec(ipt,:),opts,tspan);

    % continue using MATCONT for all npts points
    [s2,mssid,nss,hbif2] = setupMATCONT(@KotteMATCONT,@Kottecont_fluxcalc,...
                            @getKotteaxislabels,id4xeq,allpvec(ipt,:),ap,model,...
                            fluxg,npts,1500,1,hbif2); 

    % store continuation data
    cmb2siid.(['ipt' num2str(ipt)]) = s2;
    cmb2mssid.(['ipt' num2str(ipt)]) = mssid;

    ipt = ipt+4;
end
    






   


              



