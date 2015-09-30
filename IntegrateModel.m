function sol = IntegrateModel(model,pvec,mc)
%check if concentrations are initialized
if nargin < 3
    %reinitialize concentrations
    imc = zeros(model.nt_metab,1);
else
    imc = mc;
end

%initialize solver properties
[model,solverP,saveData] = imodel(model,200);

%calculate initial flux
flux = iflux(model,pvec,imc);

%integrate model
[sol,finalSS,status] = callODEsolver(model,pvec,imc,solverP);

