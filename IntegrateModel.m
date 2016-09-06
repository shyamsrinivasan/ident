function [outsol,jacobian] =...
IntegrateModel(model,ensb,mc,ess_rxn,Vup_struct,change_pos,change_neg)
% integrate iniital model to obtain initial steady state
if nargin<7
    change_neg = struct([]);
end
if nargin < 6
    change_pos = struct([]);
end
if nargin<5
    Vup_struct = ([]);
end
if nargin<4
    ess_rxn = {};
end
% check if concentrations are initialized
if nargin < 3
    error('getiest:NoA','No initial concentrations');
end
if nargin<2
    error('getiest:NoA','No parameter vector');
else
    pvec = ensb{1,2};
end

% initialize solver properties
% ecoli model
% [model,solverP,saveData] = imodel(model,1e9,ess_rxn,Vup_struct);

% toy model
[model,solverP,saveData] = imodel(model,'ode',1000,ess_rxn,Vup_struct);

% remove water and protons (held constant) from consideration in the model
% integration phase
[model,pvec,mc] = addremoveMets(model,{'h2o[c]','h2o[e]'},pvec,mc);
[newmodel,newpvec,newmc,cnstmet] = remove_eMets(model,pvec,mc,[model.Vind model.Vex],...
                           {'ac[e]','bm[e]','pep[e]'});
% metid = {'glc[e]','o2[e]','h[e]','h[c]','pi[e]','pyr[e]'}                       

% preinitialize Vind and Vex to reduce integration overhead in iflux.m                       
newmodel.Vind = addToVind(newmodel,newmodel.Vind,newmodel.rxn_add,newmodel.rxn_excep);
newmodel.rxn_excep = union(newmodel.rxn_excep,newmodel.rxns(newmodel.Vind));
newmodel.Vex = addToVind(newmodel,newmodel.Vex,[],newmodel.rxn_excep);                       

% only initialize for varmets   
nvar = length(newmodel.mets)-length(find(cnstmet));
imc = zeros(nvar,1);
imc(imc==0) = 1;

% ecoli model
% h = find(strcmpi(model.rxns,'exH'));
% model.Vuptake([h]) = [1000];

% noramlize concentration vector to intial state
Nimc = newmc(1:nvar);%imc./imc;
% Nimc(4) = 20;
PM = newmc(nvar+1:end);
Nimc(imc==0) = 0;

% intorduce perturbation in initial conditions
% met.glc_e = 10;
% Nimc(strcmpi(model.mets,'glc[e]')) = 1.1;
% if ~isempty(change_pos)
%     Nimc = changeInitialCondition(model,Nimc,change_pos);
% end
% if ~isempty(change_neg)
%     Nimc = changeInitialCondition(mdoel,Nimc,[],change_neg);
% end

newmodel.imc = imc;
newmodel.imc(newmodel.imc==0) = 1;
newmodel.PM = PM;
% calculate initial flux
flux = iflux(newmodel,newpvec,[Nimc.*imc;PM]);

% ecoli model
% dXdt = ODEmodel(0,Nimc,[],model,pvec);

% toy model
dXdt = ToyODEmodel(0,Nimc,[],newmodel,newpvec);

% get jacobian and eigen values and eigne vectors
% [J,lambda,w] = getjacobian(Nimc,newpvec,newmodel);

% integrate model
[outsol,allxeq] = callODEsolver(newmodel,newpvec,Nimc,solverP);

% introduce perturbation
Nimc = perturbEqSolution(model,allxeq.y,change_pos,change_neg);

%Perturbation to concentrations
[outsol] =...
MonteCarloPertrubationIV(model,ess_rxn,Vup_struct,ensb,allxeq.y.*imc,change_pos,[]);

% sol = MCPerturbationFlux(model,ess_rxn,Vup_struct,ensb,finalSS.y,finalSS.flux,change_pos,[]);


% initialize solver properties
[model,solverP,saveData] = imodel(model,ess_rxn,Vup_struct,1e5);

% introduce perturbation
Nimc = perturbEqSolution(model,allxeq.y,change_pos,[]);

% integrate model
[outsol,allxeq,status] = callODEsolver(model,pvec,Nimc,solverP);

% initialize solver properties
[model,solverP,saveData] = imodel(model,ess_rxn,Vup_struct,1.1e5);

% integrate model
[outsol,allxeq,status] = callODEsolver(model,pvec,Nimc,solverP,outsol);






