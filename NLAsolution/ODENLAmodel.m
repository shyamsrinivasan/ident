function sol =...
ODENLAmodel(model,ensb,mc,ess_rxn,Vup_struct,change_pos,change_neg)
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

% solution to ODE
% initialize solver properties
% toy model
[model,ODEoptions,saveData] = imodel(model,'ode',500);

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
Nimc(3) = 0.2;
Pimc = newmc(nvar+1:end);
Nimc(imc==0) = 0;
newmodel.imc = imc;
newmodel.imc(newmodel.imc==0) = 1;
newmodel.Pimc = Pimc;

% calculate initial flux
flux = iflux(newmodel,newpvec,[Nimc.*imc;Pimc]);

% toy model
dXdt = ToyODEmodel(0,Nimc,[],newmodel,newpvec);

% integrate model
% [sol,finalSS,status] = callODEsolver(newmodel,newpvec,Nimc,ODEoptions);

% solution to NLA
% initialize solver parameters for NLASS solution
% toy model
[newmodel,NLAoptions,saveData] = imodel(newmodel,'nla');

% solve NLA model
[sol,finalSS,status] = callNLAsolver(newmodel,newpvec,Nimc,NLAoptions);
