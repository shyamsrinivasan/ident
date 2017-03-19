% get volume of region of attraction based on MC simulation of initial
% values and checking final equilibrium points
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

% get unstable manifold to determine the divide between regions of
% attraction
tspanr = 0:-.1:-30;
tspanf = 0:0.1:2000;
eps = 1e-4;
[~,eig,w] = getKotteJacobian(orig_saddle,pvec,model);
[xWus,xeq] = calc1DWus(orig_saddle,w,eig,model,pvec,opts,tspanr,tspanf,eps);
% chop xWus
[~,nzid,~] = find(xWus~=0,1,'last');
relWus = real(xWus(:,1:nzid));
[x,y,z] = chopvals(relWus(1,:),relWus(2,:),relWus(3,:),[10 10 10]);

%% get region of attraction by perturbing around the steady state and 
% expanding region of perturbation
% radius where no initial point produces the 2nd steady state = rAss1 =
% [1.22 1.22 1.22]'
rndivals = get3Dsphere(xeq1',1.22,50000);
save('regionVoluemSamplevals');
% integrate
options = [];
% hfig = figure;
% [ppival,npival,allxeq,ssid,allfeq] =...
% getPequilibrium(rndivals,model,pvec,options,opts,tspanf);
% plotrndivals(rndivals,ssid,allxeq,[1 2 3],2,hfig,[]);
% plotrndivals(rndivals,ssid,allxeq,[1 2],2,hfig,[]);
% plotrndivals(rndivals,ssid,allxeq,[2 3],2,hfig,[]);
% plotrndivals(rndivals,ssid,allxeq,[2 3],2,hfig,[]);


%% get random initial values for all 3 variables and calculate steady state
% rndivals = randomivals([0 5;0 5;0 5],10000);
% % integrate
% options = [];
% hfig = figure;
% [ppival,npival,allxeq,ssid,allfeq] =...
% getPequilibrium(rndivals,model,pvec,options,opts,tspanf);
% plotrndivals(rndivals,ssid,allxeq,[1 2],2,hfig,ha)
%% get nipts random initial values and corresponding equilibrium points
% through perturbation of relWus (x,y,z)
% relWus = [x;y;z];
% newivals = perturbTrajectory(relWus);

%% integrate from newivals
% options = optimoptions('fsolve','Display','final-detailed',...
%                        'TolFun',1e-16,...
%                        'TolX',1e-12,...
%                        'MaxFunEvals',1000000,...
%                        'MaxIter',50000);
% hfig = figure;
% [allxeq,ssid,allfeq] =...
% getPequilibrium(newivals,model,pvec,options,opts,tspanf,hfig,[1 2 3]);
