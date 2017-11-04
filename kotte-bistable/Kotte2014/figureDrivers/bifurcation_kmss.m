% Figure 2B - Bifurcation diagram with stable kientic model states
% change uptake rates on line 16 (1 or 2 a.u.) to correct envelopes
runKotte

% get saddle node
[saddle,saddlepar] = getsaddlenode(data.s1,data.x1,5e-3);
pvec(ap) = saddlepar;
model.PM(ac-length(saddle)) = saddlepar;

% perturb saddle to get steady states
eps = 1e-4;
tspanf = 0:0.1:2000;
pival = saddle+eps*[1;1;1];
[~,xeq1,~,feq1] = solveODEonly(1,pival,model,pvec,opts,tspanf);
nival = saddle-eps*[1;1;1];
[~,xeq2,~,feq2] = solveODEonly(1,nival,model,pvec,opts,tspanf);
fss = [feq1 feq2];

% plot fss and saddle on bifurcation diagram - figure handle is hbif1 in
% solveEquilibriumODE
habif = get(hbif1,'Children');
set(0,'CurrentFigure',hbif1);
set(hbif1,'CurrentAxes',habif);
colorSpec = chooseColors(2,{'HotPink','Navy'});
plot(saddlepar,saddle(1),'LineStyle','none','Marker','o',...
                                            'MarkerSize',20,...
                                            'MarkerFaceColor',[.2 1 .9],...
                                            'MarkerEdgeColor',[.2 1 .9]);
plot(saddlepar,xeq1(1),'LineStyle','none','Marker','p',...
                                            'MarkerSize',20,...
                                            'MarkerFaceColor',colorSpec{2},...
                                            'MarkerEdgeColor',colorSpec{2});
plot(saddlepar,xeq2(1),'LineStyle','none','Marker','p',...
                                            'MarkerSize',20,...
                                            'MarkerFaceColor',colorSpec{1},...
                                            'MarkerEdgeColor',colorSpec{1}); 
