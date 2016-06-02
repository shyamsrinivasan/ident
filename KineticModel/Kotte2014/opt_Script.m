% optimization framework to estimate mss producing/non-producing enzyme
% parameters

% build stoichioemtrc matrices
addpath(genpath('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel'));
rxfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\Kotte2014\Kotte2014.txt';
cnfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\Kotte2014\Kotte2014C.txt';

% create model structure
[FBAmodel,parameter,variable,nrxn,nmetab] = modelgen(rxfname);

% obtain conentrations from file
[mc,FBAmodel,met] = readCNCfromFile(cnfname,FBAmodel);

% run FBA
Vup_struct.ACt2r = 1;
Vup_struct.ENZ1ex = 1;
FBAmodel = FBAfluxes(FBAmodel,'fba',{'ACt2r','ENZ1ex'},Vup_struct,...
                    [find(strcmpi(FBAmodel.rxns,'FDex'))...
                     find(strcmpi(FBAmodel.rxns,'PEPex'))]);
                 
% remove metabolites held constant from consideration in the model
% integration phase
[model,pvec,newmc,cnstmet] =...
remove_eMets(FBAmodel,parameter,mc,[FBAmodel.Vind FBAmodel.Vex],...
{'enz1[c]','enz1[e]','enz[e]','ac[e]','bm[c]','bm[e]','pep[e]'});

% only initialize for varmets   
nvar = length(model.mets)-length(find(cnstmet));
M = newmc(1:nvar);
PM = newmc(nvar+1:end);
model.PM = PM;

clear pvec
kEcat = 1;
KEacetate = 0.1;    % or 0.02
KFbpFBP = 0.1;
vFbpmax = 1;
Lfbp = 4e6;
KFbpPEP = 0.1;
vEXmax = 1;
KEXPEP = 0.3;
vemax = 1.1;        % for bifurcation analysis: 0.7:0.1:1.3
KeFBP = 0.45;       % or 0.45
ne = 2;             % or 2
acetate = 0.1;      % a.u acetate
d = 0.25;           % or 0.25 or 0.35
kPEPout = 0.2;
pvec = [KEacetate,KFbpFBP,Lfbp,KFbpPEP,...
        KEXPEP,vemax,KeFBP,ne,acetate,d,...
        kPEPout,kEcat,vFbpmax,vEXmax];
    
% calculate model jacobian
% Method 1 use ADMAT to calculate jacobians
admatfun = @(x)Kotte_givenNLAE(x,model,pvec);
% x = ones(length(M),1);
xADMATobj = deriv(M,eye(3));
xADMATres = admatfun(xADMATobj);
F = getval(xADMATres);
J = getydot(xADMATres); 

% Method 2 exact Jacobian
% Jxact = KottegivenJacobian(M,pvec,model);

% Method 3 finite difference approximation Jacobian using sparsity pattern
% get Jacobian sparsity pattern from S and SI matrices
% Jpatt = modelSparsity(model);

% obtain function handle for f(X) in Jacobian calculation
% fJrow = @(i,x)KotteSvrow(i,x,model,pvec);
% Jfd = model_Jacobian(model,[M;model.PM],fJrow,Jpatt(1:3,1:3));

% necessities for nonlinear constrained optimization
% x = [pep,fdp,enz,kEcat,vFbpmax,vEXmax];
x = [M;kEcat;vFbpmax;vEXmax]; 
pvec = [KEacetate,KFbpFBP,Lfbp,KFbpPEP,...
        KEXPEP,vemax,KeFBP,ne,acetate,d,...
        kPEPout];
    
fhandle = OptFunc(model,pvec);
% objective function  - func
obj = fhandle{1};

% gradient of objective function - grad(func)
grad = fhandle{2};

% nonlinear constraints - nlcon
nlcon = fhandle{3};

% nonlinear constraint jacobian - grad(nlcon)
nljac = fhandle{4};

% nonlinear constraint jacobian sparsity - Jpatt(nlcon)
nljacstr = fhandle{5};

% nonlinear constraint rhs - nlrhs
nlrhs = zeros(nvar,1);

% nonlinear constraint type - nle
nle = zeros(nvar,1);

% variable bounds - lb, ub
lbmet = zeros(length(M),1);
lbpar = zeros(3,1);
lbpar(lbpar==0) = 0.01;
lb = [lbmet;lbpar];
ubmet = zeros(length(M),1);
ubmet(ubmet==0) = 50;
ubpar = zeros(3,1);
ubpar(ubpar==0) = 100;
ub = [ubmet;ubpar];

% initial value - x0
x0 = [M;kEcat;vFbpmax;vEXmax]; 

% test functions obtained above from function handles
z = obj(x);
z_grad = grad(x);
A = nlcon(x);
A_jac = nljac(x);
A_jacstr = nljacstr();

% build OPTI problem object
Opt = opti('fun',obj,'grad',grad,...
           'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'nljacstr',nljacstr,...
           'bounds',lb,ub);
[x,fval,exitflag,info] = solve(Opt,x0);






