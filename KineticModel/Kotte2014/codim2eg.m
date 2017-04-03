% codim 2 bifurcation with CL_MATCONT 
% Example: catalytic oscillator
pvec = [2.5,1.92373,10,.0675,1,.1,.4];

% integrate to get initial equilibrium
odefun = @(t,x)oscillatorODE(t,x,[],pvec');
opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
tspan = 0:.1:500;
[tout,yout] = ode45(odefun,tspan,[0.00146723;0.826167;0.123119],opts);
figure
plot(tout,yout);
% [tout2,yout2] = ode45(odefun,tspan,[0.0029538;0.76211;0.16781],opts);
% figure
% plot(tout2,yout2);
% close all
ap = 2;
[data,y,p] = execMATCONT(@oscillator,[],yout(end,:)',pvec',ap,[],[],300);
fig = bifurcationPlot(data.x1,data.s1,data.f1,[4,1]);

% continue from limit points obtained on equilibrium branch
apLP = [2 7];
[LPdata_f,LPdata_b] =...
execLPcont(@oscillator,yout(end,:)',pvec',apLP,ap,data);

fig =...
bifurcationPlot(LPdata_f{1}.x1,LPdata_f{1}.s1,LPdata_f{1}.f1,[4,1],[],1,fig);
fig =...
bifurcationPlot(LPdata_b{1}.x1,LPdata_b{1}.s1,LPdata_b{1}.f1,[4,1],[],1,fig);
fig =...
bifurcationPlot(LPdata_f{2}.x1,LPdata_f{2}.s1,LPdata_f{2}.f1,[4,1],[],1,fig);
fig =...
bifurcationPlot(LPdata_b{2}.x1,LPdata_b{2}.s1,LPdata_b{2}.f1,[4,1],[],1,fig);

% calculate jacobian for points on continuation curve
% [jac,lambda,w] = getKotteJacobian(@oscillatorNLAE,yout(end,:)',p',[])
