function [sol,jacobian] =...
IntegrateModel(model,ensb,mc,ess_rxn,Vup_struct,change_pos,change_neg)
% change in initial conditions
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

%initialize solver properties
% ecoli model
% [model,solverP,saveData] = imodel(model,1e9,ess_rxn,Vup_struct);

% toy model
[model,solverP,saveData] = imodel(model,500);

% remove water and protons (held constant) from consideration in the model
% integration phase
[model,pvec,mc] = addremoveMets(model,{'h2o[c]','h2o[e]'},pvec,mc);
[newmodel,newpvec,newmc,cnstmet] = remove_eMets(model,pvec,mc,[model.Vind model.Vex],...
                           {'glc[e]','o2[e]','h[e]','h[c]','pi[e]','pyr[e]'});

% only initialize for varmets   
nvar = length(newmodel.mets)-length(find(cnstmet));
imc = zeros(nvar,1);
imc(imc==0) = 1;

% model.Vuptake = zeros(model.nt_rxn,1);


% ecoli model
% h = find(strcmpi(model.rxns,'exH'));
% model.Vuptake([h]) = [1000];

% noramlize concentration vector to intial state
Nimc = newmc(1:nvar);%imc./imc;
Pimc = newmc(nvar+1:end);
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
newmodel.Pimc = Pimc;
% calculate initial flux
flux = iflux(newmodel,newpvec,[Nimc.*imc;Pimc]);

% ecoli model
% dXdt = ODEmodel(0,Nimc,[],model,pvec);

% toy model
dXdt = ToyODEmodel(0,Nimc,[],newmodel,newpvec);

% %call to ADmat for stability/jacobian info
% [Y,Jac] = stabilityADMAT(model,pvec,Nimc.*imc);
% [e_vec,e_val] = eig(Jac);
% e_val = diag(e_val);
% jacobian = Jac;

%test reals of eigen values of jacobians
% if any(real(e_val)>0)
%     model.mets(real(e_val)>0)
% %     fprintf('%d %3.6g %d\n',find(real(e_val)>0));
% end
% Nimc_obj = deriv(Nimc,eye(model.nt_metab));
% dXdtADMAT = ODEmodelADMAT(0,Nimc_obj,[],model,pvec);
% Y = getval(dXdtADMAT);
% Jac = getydot(dXdtADMAT);

%integrate model
[sol,finalSS,status] = callODEsolver(newmodel,newpvec,Nimc,solverP);


%introduce perturbation
Nimc = perturbEqSolution(model,finalSS.y,change_pos,change_neg);

%Perturbation to concentrations
[sol] =...
MonteCarloPertrubationIV(model,ess_rxn,Vup_struct,ensb,finalSS.y.*imc,change_pos,[]);

% sol = MCPerturbationFlux(model,ess_rxn,Vup_struct,ensb,finalSS.y,finalSS.flux,change_pos,[]);


%initialize solver properties
[model,solverP,saveData] = imodel(model,ess_rxn,Vup_struct,1e5);

%introduce perturbation
Nimc = perturbEqSolution(model,finalSS.y,change_pos,[]);

%integrate model
[sol,finalSS,status] = callODEsolver(model,pvec,Nimc,solverP);

%initialize solver properties
[model,solverP,saveData] = imodel(model,ess_rxn,Vup_struct,1.1e5);

%integrate model
[sol,finalSS,status] = callODEsolver(model,pvec,Nimc,solverP,sol);

%initialize solver properties
% [model,solverP,saveData] = imodel(model,ess_rxn,Vup_struct,5e3);

%integrate model
% [sol,finalSS,status] = callODEsolver(model,pvec,Nimc,solverP);

%initialize solver properties
[model,solverP,saveData] = imodel(model,ess_rxn,Vup_struct,3e5);

%integrate model
[sol,finalSS,status] = callODEsolver(model,pvec,Nimc,solverP,sol);

%initialize solver properties
% [model,solverP,saveData] = imodel(model,ess_rxn,Vup_struct,5e4);

%integrate model
% [sol,finalSS,status] = callODEsolver(model,pvec,Nimc,solverP);

%initialize solver properties
[model,solverP,saveData] = imodel(model,ess_rxn,Vup_struct,4e5);

%integrate model
[sol,finalSS,status] = callODEsolver(model,pvec,Nimc,solverP,sol);

%initialize solver properties
[model,solverP,saveData] = imodel(model,ess_rxn,Vup_struct,5e5);

%integrate model
[sol,finalSS,status] = callODEsolver(model,pvec,Nimc,solverP,sol);

%initialize solver properties
[model,solverP,saveData] = imodel(model,ess_rxn,Vup_struct,7e5);

%integrate model
[sol,finalSS,status] = callODEsolver(model,pvec,Nimc,solverP,sol);

%initialize solver properties
[model,solverP,saveData] = imodel(model,ess_rxn,Vup_struct,9e5);

%integrate model
[sol,finalSS,status] = callODEsolver(model,pvec,Nimc,solverP,sol);

%initialize solver properties
[model,solverP,saveData] = imodel(model,ess_rxn,Vup_struct,1e6);

%integrate model
[sol,finalSS,status] = callODEsolver(model,pvec,Nimc,solverP,sol);




