function [sol,jacobian] =...
IntgegrateNLASSmodel(model,ensb,mc,ess_rxn,Vup_struct,change_pos,change_neg)
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
if nargin < 3
    error('getiest:NoA','No initial concentrations');
end
if nargin<2
    error('getiest:NoA','No parameter vector');
else
    pvec = ensb{1,2};
end

% initialize solver parameters for NLASS solution
% toy model
[model,solverP,saveData] = imodel(model,'nla');

% adjust model for solution
% remove water and protons (held constant) from consideration in the model
% integration phase
[model,pvec,mc] = addremoveMets(model,{'h2o[c]','h2o[e]'},pvec,mc);
[newmodel,newpvec,newmc,cnstmet] = remove_eMets(model,pvec,mc,[model.Vind model.Vex],...
                           {'glc[e]','o2[e]','h[e]','h[c]','pi[e]','pyr[e]'});

% only initialize for varmets   
nvar = length(newmodel.mets)-length(find(cnstmet));
imc = zeros(nvar,1);
imc(imc==0) = 1;    

% noramlize concentration vector to intial state
Nimc = newmc(1:nvar);%imc./imc;
Pimc = newmc(nvar+1:end);
Nimc(imc==0) = 0;
newmodel.imc = imc;
newmodel.imc(newmodel.imc==0) = 1;
newmodel.Pimc = Pimc;

% calculate initial flux
flux = iflux(newmodel,newpvec,[Nimc.*imc;Pimc]);

% calculate rhs of 0 = dx/dt = f(x)
dXdt = ToyODEmodel(0,Nimc,[],newmodel,newpvec);

% solve NLA model
[sol,finalSS,status] = callNLAsolver(newmodel,newpvec,Nimc,solverP);



