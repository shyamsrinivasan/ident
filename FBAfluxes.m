function model = FBAfluxes(model)

if ~isfield(model,'Vuptake')
    Vuptake = zeros(model.nt_rxn,1);
    Vuptake(strcmpi(model.rxns,'exGLC')) = 20;%mmol/gDCW.h
    Vuptake(strcmpi(model.rxns,'exO2')) = 10;%mmole/gDCW.h
%         Vuptake(strcmpi(model.rxns,'exH2O')) = 40;
%         Vuptake(strcmpi(model.rxns,'exPI')) = 50;
    model.Vuptake = Vuptake;
end
    
%calculate FBA fluxes and growth rate
[~,vLP,flag] = estimateLPgrowth(model);
if flag
    %steady state FBA fluxes for direction
    model.Vss = vLP;
else
    fprintf('FBA fluxes unassigned\n');
end

    