% find points on the stability boundary for Kotte model
% Chiang, Hirsch and Wu, 1988 (refered from Khalil, Nonlinear Systems)
% step 1. find all equilibrium points
runKotte

% get saddle node to get both stable steady states and get a bistable
% acetate concentration for perturbation
[orig_saddle,orig_saddlepar,saddleid] = getsaddlenode(data.s1,data.x1,5e-3);
pvec(ap) = orig_saddlepar;
model.PM(ac-length(orig_saddle)) = orig_saddlepar;
saddleflux = Kotte_givenFlux([orig_saddle;model.PM],pvec,model);
% get eigen values and eigen vectors of the saddle node
[~,alleig,alleigw] = KotteStabilityInfo(orig_saddle',model,pvec);

% perturb saddle to get steady states
eps = 1e-4;
tspanf = 0:0.1:2000;
pival = orig_saddle+eps*[1;1;1];
[~,xeq1,~,feq1] = solveODEonly(1,pival,model,pvec,opts,tspanf);
nival = orig_saddle-eps*[1;1;1];
[~,xeq2,~,feq2] = solveODEonly(1,nival,model,pvec,opts,tspanf);
xss = [xeq1 xeq2];

%% pick points on the boundary of region of attraction - unstable eigen vectors
figure
hold on
itermax = 100;
eps = 1e-3;
alpha = 0.8;
[ip_interUS,status] =...
getUSintsecpt(orig_saddle,alleigw(:,3),model,pvec,opts,-1,eps,alpha,itermax);
for ip = 1:size(ip_interUS,2)
    xfwd = solveODEonly(1,ip_interUS(:,ip),model,pvec,opts,[0,2000]);
    xrev = solveODEonly(1,ip_interUS(:,ip),model,pvec,opts,[0,-40]);
%     plot(ip_interUS(1,ip),ip_interUS(2,ip),'Marker','.','MarkerSize',10,'Color','b');
%     plot(xfwd(1,:),xfwd(2,:),'k');
%     plot(xrev(1,:),xrev(2,:),'r');
    plot3(ip_interUS(1,ip),ip_interUS(2,ip),ip_interUS(3,ip),'Marker','.','MarkerSize',10,'Color','b');
    plot3(xfwd(1,:),xfwd(2,:),xfwd(3,:),'k');
    drawnow
    plot3(xrev(1,:),xrev(2,:),xrev(3,:),'r');
    drawnow
end

%% pick points on the boundary of region of attraction - stable eigen vectors
% figure
% hold on
itermax = 5;
eps = 1e-3;
alpha = 0.8;
for iw = 1:2
    [ip_interS,status] =...
    getSintsecpt(orig_saddle,alleigw(:,iw),model,pvec,opts,[1 2],eps,alpha,itermax);
    for ip = 1:size(ip_interUS,2)
        xfwd = solveODEonly(1,ip_interS(:,ip),model,pvec,opts,[0,2000]);
        xrev = solveODEonly(1,ip_interS(:,ip),model,pvec,opts,[0,-40]);
    %     plot(ip_interS(1,ip),ip_interS(2,ip),'Marker','.','MarkerSize',10,'Color','b');
    %     plot(xfwd(1,:),xfwd(2,:),'k');
    %     plot(xrev(1,:),xrev(2,:),'r');
        plot3(ip_interS(1,ip),ip_interS(2,ip),ip_interS(3,ip),'Marker','.','MarkerSize',10,'Color','b');
        plot3(xfwd(1,:),xfwd(2,:),xfwd(3,:),'k');
        drawnow
        plot3(xrev(1,:),xrev(2,:),xrev(3,:),'r');
        drawnow
    end
end