% algorithm for calculating the approximate region of attraction of
% dyamical systems
% Chiang, Hirsch and Wu, 1988 (refered from Khalil, Nonlinear Systems)

% step 1. find all equilibrium points
runKotte

% get saddle node to get both stable steady states and get a bistable
% acetate concentration for perturbation
[orig_saddle,orig_saddlepar,saddleid] = getsaddlenode(data.s1,data.x1,5e-3);
pvec(ap) = orig_saddlepar;
model.PM(ac-length(orig_saddle)) = orig_saddlepar;
saddleflux = Kotte_givenFlux([orig_saddle;model.PM],pvec,model);

% perturb saddle to get steady states
eps = 1e-4;
tspanf = 0:0.1:2000;
pival = orig_saddle+eps*[1;1;1];
[~,xeq1,~,feq1] = solveODEonly(1,pival,model,pvec,opts,tspanf);
nival = orig_saddle-eps*[1;1;1];
[~,xeq2,~,feq2] = solveODEonly(1,nival,model,pvec,opts,tspanf);
xss = [xeq1 xeq2];
                                                                 
% all equilibrium points & corresponding eigen values
% step 2: get unit length eigen vectors of unstable eq. pt.
alleqpts = [xss orig_saddle]';
[~,alleig,alleigw] = KotteStabilityInfo(alleqpts,model,pvec);
unsteig = alleig(:,3);
unsteigw = alleigw((3-1)*nvar+1:3*nvar,:);
% Eus = unsteigw(:,real(unsteig)>0);

calcBasinBoundary(xss,unsteig,unsteigw,model,pvec,opts,0:-20,tspanf,5e-3)





