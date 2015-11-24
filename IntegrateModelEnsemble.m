function [sol_all,Jac_all] = IntegrateModelEnsemble(model,ess_rxn,Vup_struct,ensb,mc,change_pos,change_neg)
%same concentration vector mc is used for all models in ensb

%change in initial conditions
if nargin<7
    change_neg = struct([]);
end
if nargin < 6
    change_pos = struct([]);
end
%check if concentrations are initialized
if nargin < 5
    %reinitialize concentrations
    imc = zeros(model.nt_metab,1);
else
    imc = mc;
end

if nargin<4
    error('getiest:NoA','No parameter vector');       
end

if nargin<3
    Vup_struct([]);
end
if nargin<2
    ess_rxn = {};
end

%initialize required variables
nmodels = size(ensb,1);
ini_flux = cell(nmodels,1);
ini_dXdt = cell(nmodels,1);
Jac_all = cell(nmodels,1);
e_val = cell(nmodels,1);

sol_all = cell(nmodels,1);
finalSS_all = cell(nmodels,1);

%initialize solver properties
%once for all models
[model,solverP,saveData] = imodel(model,ess_rxn,Vup_struct,1.2);

%noramlize concentration vector to intial state
Nimc = imc./imc;
Nimc(imc==0) = 0;

%intorduce perturbation in initial conditions
% met.glc_e = 10;
% Nimc(strcmpi(model.mets,'glc[e]')) = 1.1;
if ~isempty(change_pos)
    Nimc = changeInitialCondition(model,Nimc,change_pos);
end
if ~isempty(change_neg)
    Nimc = changeInitialCondition(mdoel,Nimc,[],change_neg);
end

model.imc = imc;
model.imc(model.imc==0) = 1;

for im = 1:nmodels
    if ensb{im,2}.feasible
        pvec = ensb{im,2};
    else
        fprintf('Model %d in ensemble is infeasible\n',im);
        fprintf('Continuing to next model: Model %d\n',im);
        continue;
    end
    
    %calculate initial flux
    flux = iflux(model,pvec,Nimc.*imc);
    dXdt = ODEmodel(0,Nimc,[],model,pvec);
    ini_flux{im} = flux;
    ini_dXdt{im} = dXdt;
    
    % %call to ADmat for stability/jacobian info
    [Y,Jac] = stabilityADMAT(model,pvec,Nimc);
    e_val{im} = eig(Jac);
    Jac_all{im} = Jac;
    
    %test reals of eigen values of jacobians
    if any(real(e_val{im})>0)
        model.mets(real(e_val{im})>0)    
    end
    
    %integrate model
    [sol,finalSS,status] = callODEsolver(model,pvec,Nimc,solverP);
    sol_all{im} = sol;
    finalSS_all{im} = finalSS;    
end