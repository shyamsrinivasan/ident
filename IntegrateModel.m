function sol = IntegrateModel(model,ensb,mc)
%check if concentrations are initialized
if nargin < 3
    %reinitialize concentrations
    imc = zeros(model.nt_metab,1);
%     imc(strcmpi(model.mets,'pep[c]')) = 0.002;
%     imc(strcmpi(model.mets,'pep[c]')) = 1e-5;
else
    imc = mc;
%     imc(strcmpi(model.mets,'h[c]')) = 100;
%     imc(strcmpi(model.mets,'q8h2[c]')) = 100;
%     imc(strcmpi(model.mets,'pi[c]')) = 100;
end

if nargin<2
    error('getiest:NoA','No parameter vector');
else
    pvec = ensb{1,2};
end
    
%initialize solver properties
[model,solverP,saveData] = imodel(model,100);

% model.Vuptake = zeros(model.nt_rxn,1);
% h2o = find(strcmpi(model.rxns,'exH2O'));
% pi =  find(strcmpi(model.rxns,'exPI'));
h = find(strcmpi(model.rxns,'exH'));
% 
model.Vuptake([h]) = [1000];

%noramlize concentration vector to intial state
Nimc = imc./imc;
Nimc(imc==0) = 0;
model.imc = imc;
model.imc(model.imc==0) = 1;
%calculate initial flux
flux = iflux(model,pvec,Nimc.*imc);
dXdt = ODEmodel(0,Nimc,[],model,pvec);

%integrate model
[sol,finalSS,status] = callODEsolver(model,pvec,Nimc,solverP);

