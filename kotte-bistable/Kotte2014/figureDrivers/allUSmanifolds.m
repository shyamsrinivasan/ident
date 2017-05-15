runKotte 

% get saddle node
[saddle,saddlepar] = getsaddlenode(data.s1,data.x1,5e-3);
% get all parameter (saddles) values between the 2 limit points
allsaddles = extractsaddles(data.s1,data.x1);
nsadpts = size(allsaddles,2);

% get eigval and eig vec for all saddle points
alleig = zeros(nvar,nsadpts);
allw = zeros(nvar,nvar,nsadpts);
tspanr = 0:-.1:-30;
tspanf = 0:0.1:2000;
eps = 1e-4;
allWus = zeros(nvar,length(tspanr),nsadpts);

for ipts = 1:nsadpts
    pvec(ap) = allsaddles(end,ipts);
    model.PM(ac-length(saddle)) = allsaddles(end,ipts);
    
    % get eig val and eig vector
    [~,alleig(:,ipts),w] = getKotteJacobian(allsaddles(1:3,ipts),pvec,model);
    allw(:,:,ipts) = w;
    
    % calculate unstable manifold for all saddle points
    [xWus,xeq] = calc1DWus(allsaddles(1:3,ipts),w,alleig(:,ipts),model,pvec,opts,tspanr,tspanf,eps);
    
    % collect reverse trajectory 
    allWus(:,:,ipts) = xWus;
end

%% unstable manifold in 3D
% figure
allx = [];ally = [];allz = [];
for jpt = 1:nsadpts    
    [~,nzid,~] = find(allWus(:,:,jpt)~=0,1,'last');
    relWus = allWus(:,1:nzid,jpt);
    kid = relWus(1,:)<0|relWus(1,:)>5|relWus(2,:)<0|relWus(2,:)>5|relWus(3,:)<0|relWus(3,:)>5;
    relWus(:,kid) = [];
    allx = [allx relWus(1,:)];
    ally = [ally relWus(2,:)];
    allz = [allz relWus(3,:)];
%     hold on
%     plot3(relWus(1,:),relWus(2,:),relWus(3,:),'LineStyle','none','Marker','.','Color','r');    
end
hfig = Manifold2DPlot(allx,ally,allz);


% perturbManifolds(xmafold(1:3,:)',saddle,model,pvec,opts,tspanf,...
%                         hfig,[1 2 3],10); 