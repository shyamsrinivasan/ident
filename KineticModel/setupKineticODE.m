function [newmodel,newpvec,Nimc,solverP,flux,dXdt] =...
        setupKineticODE(model,ensb,mc,ess_rxn,Vup_struct,tmax)
% setup kinetic model for integration
if nargin<6
    tmax = 500;
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
[model,solverP,saveData] = imodel(model,'ode',tmax,ess_rxn,Vup_struct);

% remove water and protons (held constant) from consideration in the model
% integration phase
[model,pvec,mc] = addremoveMets(model,{'h2o[c]','h2o[e]'},pvec,mc);
[newmodel,newpvec,newmc,cnstmet] = remove_eMets(model,pvec,mc,[model.Vind model.Vex],...
                           {'ac[e]','bm[e]','pep[e]'});
                       
% preinitialize Vind and Vex to reduce integration overhead in iflux.m                       
newmodel.Vind = addToVind(newmodel,newmodel.Vind,newmodel.rxn_add,newmodel.rxn_excep);
newmodel.rxn_excep = union(newmodel.rxn_excep,newmodel.rxns(newmodel.Vind));
newmodel.Vex = addToVind(newmodel,newmodel.Vex,[],newmodel.rxn_excep);  

% only initialize for varmets   
nvar = length(newmodel.mets)-length(find(cnstmet));
imc = zeros(nvar,1);
imc(imc==0) = 1;

% separate variables from constants and assign constants to a field in
% newmodel
Nimc = newmc(1:nvar);
PM = newmc(nvar+1:end);
Nimc(imc==0) = 0;
newmodel.imc = imc;
newmodel.imc(newmodel.imc==0) = 1;
newmodel.PM = PM;

% system check
% calculate initial flux
flux = iflux(newmodel,newpvec,[Nimc.*imc;PM]);
% toy model
dXdt = ToyODEmodel(0,Nimc,[],newmodel,newpvec);
% get jacobian and eigen values and eigne vectors
% [J,lambda,w] = getjacobian(Nimc,newpvec,newmodel);
