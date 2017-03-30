% get new unstable manifolds for different enzyme expression levels
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

% run enzyme perturbations
% sample parameters indicated by indices in idp
cmb = [.05 1 1;1 .05 1;1 1 .05;.05 .05 .05;...
       .125 1 1;1 .125 1;1 1 .125;.125 .125 .125;...
       .25 1 1;1 .25 1;1 1 .25;.25 .25 .25;...
       .5 1 1;1 .5 1;1 1 .5;.5 .5 .5;...
       2 1 1;1 2 1;1 1 2;2 2 2;...
       4 1 1;1 4 1;1 1 4;4 4 4];

idp = [12 13 14];
type = 'together';
npts = size(cmb,1);

if strcmpi(type,'together')
    alliidpvec = zeros(npts,length(pvec),size(idp,1));    
else    
    alliidpvec = zeros(npts,length(pvec),length(idp));    
end

% set acetate conentration
pvec = [K1ac,K3fdp,L3,K3pep,...
        K2pep,vemax,KeFBP,ne,acetate,d,...
        k4cat,k1cat,v3max,v2max];
pvec(ap) = 0.01;
model.PM(ac-length(orig_saddle)) = 0.01;

% set parameters from cmb at idp position(s)
allpvec = repmat(pvec,npts,1);
allpvec(:,idp) = cmb;

% find equilibrium point for lowest acetate 
[~,allxeqlac] = solveODEonly(npts,M,model,allpvec,opts,tspan);

% and continue on acetate
[s,mssid,nss] = setupMATCONT(@KotteMATCONT,@Kottecont_fluxcalc,allxeqlac,...
                            allpvec,ap,model,fluxg,npts,1500);

%% extract all unstable points for each enzyme perturbation that is bistable
% figure
% hold on
for ipt = 1:npts
    if ismember(ipt,mssid)
        allsaddles = extractsaddles(s.(['pt' num2str(ipt)]).s1,s.(['pt' num2str(ipt)]).x1);       
        nsadpts = size(allsaddles,2);
%         plot3(allsaddles(1,:),allsaddles(2,:),allsaddles(3,:),'Color','r');
    end
end

% Wus for every point on the unstable bifurcation curve
allWus2 = cell(length(mssid),3);
tspanr = 0:-.1:-30;
tspanf = 0:0.1:2000;
eps = 1e-4;
iss = 1;
figure
hold on
for ipt = 1:npts
    if ismember(ipt,mssid)
        allx = [];ally = [];allz = [];        
        ispvec = allpvec(ipt,:);
        allsaddles = extractsaddles(s.(['pt' num2str(ipt)]).s1,s.(['pt' num2str(ipt)]).x1);       
        nsadpts = size(allsaddles,2);
        alleig = zeros(nvar,nsadpts);
        allw = zeros(nvar,nvar,nsadpts);        
        for isad = 1:nsadpts
            ispvec(ap) = allsaddles(end,isad);
            model.PM(ac-length(orig_saddle)) = allsaddles(end,isad);
            % get eig val and eig vector
            [~,alleig(:,isad),w] = getKotteJacobian(allsaddles(1:3,isad),ispvec,model);
            allw(:,:,isad) = w;
            % calculate unstable manifold for all saddle points
            [xWus,xeq] =...
            calc1DWus(allsaddles(1:3,isad),w,alleig(:,isad),model,ispvec,opts,tspanr,tspanf,eps);
            % collect reverse trajectory 
            if ~isempty(xWus)                
                % chop and plot all points in allWus2
                [~,nzid] = find(xWus~=0,1,'last');
                xWus = real(xWus(:,1:nzid));
                [x,y,z] = chopvals(xWus(1,:),xWus(2,:),xWus(3,:),[20 20 20]);
                allx = [allx x];ally = [ally y];allz = [allz z];
            end
        end
        allWus2{iss,1} = allx;
        allWus2{iss,2} = ally;
        allWus2{iss,3} = allz;
        iss = iss+1;
    end
end

%% get boundaries of acetate bistability
allEbsval = allpvec(mssid,idp);
acbounds = zeros(2,length(mssid)); % [min;max];
% xbounds = zeros(nvar,2*length(mssid));
mssipt = 1;
for ipt = 1:npts
    if ismember(ipt,mssid)
        index = cat(1,s.(['pt' num2str(ipt)]).s1.index);
        x1 = s.(['pt' num2str(ipt)]).x1;
%         xcont = x1(1:nvar,index);
        pcont = x1(nvar+1:end,index);
%         xbounds(:,2*mssipt-1:2*mssipt) = xcont(:,2:end-1);
        acbounds(1,mssipt) = min(pcont(:,2:end-1));
        acbounds(2,mssipt) = max(pcont(:,2:end-1));
        mssipt = mssipt+1;
    end
end

%% get a saddle node for each of the bistable cases and plot Wus
% alleig = zeros(nvar,length(mssid));
% allw = zeros(nvar,nvar,length(mssid));
% allsaddles = zeros(nvar,length(mssid));
% allxeq = zeros(nvar,length(mssid));
% issid = 1;
% tspanr = 0:-.1:-30;
% tspanf = 0:0.1:2000;
% eps = 1e-4;
% eps1 = 1e-4;
% allWus = zeros(nvar,length(tspanr),length(mssid));
% % pick a saddle node    
% [saddle,saddlepar,status] = geteqcontpt(s,orig_saddlepar,eps);
% 
% figure
% hold on
% for iss = 1:npts
%     if ismember(iss,mssid)
%         ispvec = allpvec(iss,:);
%                 
%         % get eig val and eig vector
%         ispvec(ap) = saddlepar(iss);
%         model.PM(ac-nvar) = saddlepar(iss);
%         [~,alleig(:,issid),w] = getKotteJacobian(saddle(:,iss),ispvec,model);
%         allw(:,:,issid) = w;
%         
%         % get Wus at saddlepar        
%         [xWus,xeq] =...
%         calc1DWus(saddle(:,iss),w,alleig(:,issid),model,ispvec,opts,tspanr,tspanf,eps1);
%         
%         % collect reverse trajectory 
%         allWus(:,:,issid) = xWus;
%         allxeq(:,issid) = xeq(:,1);      
%         
%         issid = issid+1;
%     end
% end
% 
% % unstable manifolds after chopping
% relallWus = real(allWus);
% figure
% hold on
% for jpt = 1:length(mssid)
%     [~,nzid,~] = find(relallWus(:,:,jpt)~=0,1,'last');
%     relWus = relallWus(:,1:nzid,jpt);
%     kid = relWus(1,:)<0|relWus(1,:)>20|relWus(2,:)<0|relWus(2,:)>20|relWus(3,:)<0|relWus(3,:)>20;
%     relWus(:,kid) = [];
% %     allx = [allx relWus(1,:)];
% %     ally = [ally relWus(2,:)];
% %     allz = [allz relWus(3,:)];
%     % plot trajectory
%     plot3(relWus(1,:),relWus(2,:),relWus(3,:),'LineWidth',2);
%     plot3(saddle(1,mssid(jpt)),saddle(2,mssid(jpt)),saddle(3,mssid(jpt)),'Marker','.','Color','r','MarkerSize',16);
%     drawnow
% end



%%
% colorSpec = chooseColors(4,{'Navy','HotPink','Red','Orange'});
% for iid = 1:1 % length(idp)
%     fprintf('Parameter Combination #%d\n',iid); 
%     % find equilibrium solution from different ss at the original saddle
%     % node
%     alliidpvec(:,:,iid) = allpvec;    
%     hf1 = [];
%     ha1 = [];
%     hf2 = [];
%     ha2 = [];
% 
%     ivalpts = zeros(2*nvar,npts);
%     xeqpts = zeros(2*nvar,npts);
%     eqid = zeros(2,npts);
%     ivalid = zeros(2,npts);
%     % perturbation for all points
%     for ipt = 1:npts
%         pvec = alliidpvec(ipt,:,iid);                
%         % perturbations from ss (xeq1 and xeq2)
%         [ivalpts,ivalid,xeqpts,eqid,hf1,ha1] = ParameterPerturbations(model,pvec,...
%             xss,ivalpts,ivalid,xeqptseqid,ipt,tspanf,colorSpec,opts,hf1,ha1);                          
%     end    
% end