function [Y,Jac] = stabilityADMAT(model,ensb,mc)
if nargin < 3
    imc = zeros(model.nt_metab,1);
else
    imc = mc;
end

if nargin<2
    error('getiest:NoA','No parameter vector');
else
    pvec = ensb{1,2};
end

Nimc = imc;%./imc;
Nimc(imc==0) = 0;

model.imc = ones(model.nt_metab,1);%imc;
model.imc(model.imc==0) = 1;

%call to check actual dXdt
%calculate initial flux
flux_ = iflux(model,pvec,Nimc.*imc);
dXdt_ = ODEmodel(0,Nimc,[],model,pvec);

% %call to admat for stability/jacobian info
Nimc_obj = deriv(Nimc,eye(model.nt_metab));
dXdtADMAT = ODEmodelADMAT(0,Nimc_obj,[],model,pvec);
Y = getval(dXdtADMAT);
Jac = getydot(dXdtADMAT);

%change model - rem0ove reactions
% [model,pvec] = removeRxns(model,pvec,model.Vex);
% [model,pvec] = removeRxns(model,pvec,model.bmrxn);
% [model,pvec] = removeRxns(model,pvec,find(model.Vss==0));
% model.S = model.S(:,[1:7]);
% model.Vex = setdiff(model.Vex,1:7);

% flux = iflux(model,pvec,Nimc.*imc);
% dXdt = ODEmodel(0,Nimc,[],model,pvec);

% model.S = model.S(:,[1:7]);
% model.Vex = setdiff(model.Vex,1:7);
% %call to ADmat for stability/jacobian info
% Nimc_obj = deriv(Nimc,eye(model.nt_metab));
% dXdtADMAT = ODEmodelADMAT(0,Nimc_obj,[],model,pvec);
% Y = getval(dXdtADMAT);
% Jac = getydot(dXdtADMAT);