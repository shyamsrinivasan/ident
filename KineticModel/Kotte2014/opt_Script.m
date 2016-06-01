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
kPEPout = 2.0004e-4;
pvec = [kEcat,KEacetate,...
        KFbpFBP,vFbpmax,Lfbp,KFbpPEP,...
        vEXmax,KEXPEP,...
        vemax,KeFBP,ne,acetate,d,kPEPout];
    
% calculate model jacobian
% Method 1 use ADMAT to calculate jacobians
admatfun = @(x)Kotte_givenNLAE(x,model,pvec);
% x = ones(length(M),1);
xADMATobj = deriv(M,eye(3));
xADMATres = admatfun(xADMATobj);
F = getval(xADMATres);
J = getydot(xADMATres); 

% Method 2 exact Jacobian
Jxact = KottegivenJacobian(M,pvec,model);

% Method 3 finite difference approximation Jacobian
% get Jacobian sparsity pattern from S and SI matrices
Jpatt = modelSparsity(model);

% obtain function handle for f(X) in Jacobian calculation
fJrow = @(i,x)KotteSvrow(i,x,model,pvec);

Jfd = model_Jacobian(model,[M;model.PM],fJrow,Jpatt(1:3,1:3));

% necessities for nonlinear constrained optimization
fhandle = OptFunc(x,model);
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
lb
% initial value - x0

cplexlp
Opt = opti('fun',fun,'bounds',lb,ub,'x0',x0);

Opt = opti('f',f,'ineq',A,b,'bounds',lb,ub);
% f = c;
% A = [a11 a12;a21 a22];
% b = [0;0];
% lb = [0;0];
% ub = [0;0];
[x,fval,exitflag,info] = solve(Opt);

% fun - nonlinear objective function
Opt = opti('f',f,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'nljacstr',nljacstr);

[x,fval,exitflag,info] = solve(Opt,x0);

% use ADMAT to calculate jacobians
fun = ADfun('',length(M));
y = deriv(M,diag(ones(3,1)));


