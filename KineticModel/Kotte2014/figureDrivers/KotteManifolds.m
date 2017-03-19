% get different stable maniforlds from Kotte model using information from
% Lyons et al., 2014, Int. J. Bifurcation Chaos
% plot both unstable and stable manifolds for a given saddle point

% generate equilibrium solution and model for Kotte model
runKotte

% get saddle node
[saddle,saddlepar] = getsaddlenode(data.s1,data.x1,5e-3);
pvec(ap) = saddlepar;
model.PM(ac-length(saddle)) = saddlepar;

% calculate eig val for all points on continuation curve data.x1
alleig = zeros(3,size(data.x1,2));
for ipts = 1:size(data.x1,2)
    pvec(ap) = data.x1(end,ipts);
    model.PM(ac-length(saddle)) = data.x1(end,ipts);
    [~,alleig(:,ipts)] = KotteStabilityInfo(data.x1(1:3,ipts)',model,pvec);      
end

% perturb saddle to get steady states
pvec(ap) = saddlepar;
model.PM(ac-length(saddle)) = saddlepar;
eps = 1e-4;
tspanf = 0:0.1:2000;
pival = saddle+eps*[1;1;1];
[~,xeq1] = solveODEonly(1,pival,model,pvec,opts,tspanf);

nival = saddle-eps*[1;1;1];
[~,xeq2] = solveODEonly(1,nival,model,pvec,opts,tspanf);

[~,eigval,w] = getKotteJacobian(saddle,pvec,model);

% all manifolds in 2D
% tspanr = [0,-30;0 -12];
% NumericalSeparatrix(model,pvec,opts,ap,data.s1,data.x1,[xeq1 xeq2],'all',tspanr,2,5e-3);

% unstable manifold in 3D
tspanr = [0,-30]; % 0 -12];
hfig =...
NumericalSeparatrix(model,pvec,opts,ap,data.s1,data.x1,[xeq1 xeq2],'unstable',tspanr,3,5e-3);

% Lyons et al., 2014 code
% mypath = 'C:\Users\shyam\Documents\MATLAB\CCFM_manifolds';
% addpath(strcat(mypath,'\CCFM_manifolds\functions\'));

% set viewpoint for 3D plots
azimuth = 52; elevation = 10;

% color_hash for plotting
color_hash = {'k.','go','ko'};

eqpts = [saddle'; % 2D stable 
       xeq1';   % spiral sink (3D stable)
       xeq2'];  % sink (3D stable)

npoints = 500;
range = [saddlepar 2];
contdir = 1;
options = optimoptions('fsolve','Display','off','TolFun',1e-10,'TolX',1e-10);
% fixed_points = kotte_branches(npoints,range,contdir,eqpts,model,pvec,ap,options);   
[type,alleig] = KotteStabilityInfo(eqpts,model,pvec);      

% get points on 2D stable invariant manifold surface
tspanr = 0:-.1:-25;
[xchop,ychop,zchop] = calc2DWs(saddle,saddlepar,ap,model,pvec,tspanr,opts);


%% Separate out segments of manifold surface
% working with xchop
% backup 
x1 = real(xchop);y1 = real(ychop);z1 = real(zchop);
[xs1,ys1,zs1] = removeredundantpoints(x1,y1,z1,0.0001);
hfig = Manifold2DPlot(xs1,ys1,zs1,hfig);

% id2 = y1<1.1885;
% id1 = setdiff(1:length(x1),find(id2));
% 
% xs1 = x1(id1);
% ys1 = y1(id1);
% zs1 = z1(id1);
% x1(id1) = [];y1(id1) = [];z1(id1) = [];
% 
% % x1(id1) = [];y1(id1) = [];z1(id1) = [];
% % [x1,y1,z1] = removeredundantpoints(real(x1),real(y1),real(z1),0.01);
% 
% id3 = x1<1.325&y1>0.485&z1<2.6065;
% xs2 = x1(id3);ys2 = y1(id3);zs2 = z1(id3);
% x1(id3) = [];y1(id3) = [];z1(id3) = [];
% % 
% % % small segmet in xs2
% id3 = xs2>.2491 & xs2<.4319 & ys2>1.0905 & ys2<1.195 & zs2>.5164 & zs2<.5796;
% xs1 = [xs1 xs2(id3)];ys1 = [ys1 ys2(id3)];zs1 = [zs1 zs2(id3)];
% [xs1,ys1,zs1] = removeredundantpoints(xs1,ys1,zs1,0.02);
% xs2(id3) = [];ys2(id3) = [];zs2(id3) = [];
% 
% id4 = x1<0.687 & z1<2.4;
% xs3 = x1(id4);ys3 = y1(id4);zs3 = z1(id4);
% x1(id4) = [];y1(id4) = [];z1(id4) = [];
% % divide xs3 into 2 parts
% id5 = ys3>0.06834;
% id6 = setdiff(1:length(ys3),find(id5));
% xs31 = xs3(id5);ys31 = ys3(id5);zs31 = zs3(id5);
% xs32 = xs3(id6);ys32 = ys3(id6);zs32 = zs3(id6);
% 
% hfig = Manifold2DPlot(xs1,ys1,zs1,hfig);
% Manifold2DPlot(xs2,ys2,zs2,hfig);
% Manifold2DPlot(xs31,ys31,zs31,hfig);
% Manifold2DPlot(xs32,ys32,zs32,hfig);
% Manifold2DPlot(x1,y1,z1,hfig);

gcf
axis([0 2.5 0 1.6 0 1.6]);
view([116 22]);
fname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\Results\KotteModel\manifolds\manifoldFig';
print('-depsc','-painters','-loose',fname)

%% perturb 2D stable manifold
pvec(ap) = saddlepar;
model.PM(ac-length(saddle)) = saddlepar;
perturbManifolds([xs1' ys1' zs1'],saddle,model,pvec,opts,tspanf,hfig,3);

% plot3(saddle(1),saddle(2),saddle(3),'Color','r','Marker','.','MarkerSize',30);

% 
% % foil on surface xs2
% id5 = xs2<1.325 & ys2>0.7312 & ys2<1.189 & zs2>0.6305 & zs2<2.606;
% xs3 = xs2(id5);ys3 = ys2(id5);zs3 = zs2(id5);
% xs2(id5) = [];ys2(id5) = [];zs2(id5) = [];

% id6 = xs2>0.0004287 & xs2<0.3842 & ys2>0.7312 & ys2<0.9922 & zs2>0.5523 & zs2<0.61375;
% xs3 = [xs3 xs2(id6)];ys3 = [ys3 ys2(id6)];zs3 = [zs3 zs2(id6)];
% xs2(id6) = [];ys2(id6) = [];zs2(id6) = [];                    
                       