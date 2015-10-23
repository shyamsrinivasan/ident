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
[model,solverP,saveData] = imodel(model,10);

% model.Vuptake = zeros(model.nt_rxn,1);
% h2o = find(strcmpi(model.rxns,'exH2O'));
% pi =  find(strcmpi(model.rxns,'exPI'));
h = find(strcmpi(model.rxns,'exH'));
% 
model.Vuptake([h]) = [1000];

%calculate initial flux
flux = iflux(model,pvec,imc);


%integrate model
[sol,finalSS,status] = callODEsolver(model,pvec,imc,solverP);

