function model = FBAfluxes(model,option)

if ~isfield(model,'Vuptake')
    Vuptake = zeros(model.nt_rxn,1);
    Vuptake(strcmpi(model.rxns,'exGLC')) = 20;%mmol/gDCW.h
    Vuptake(strcmpi(model.rxns,'exO2')) = 100;%mmole/gDCW.h
%         Vuptake(strcmpi(model.rxns,'exH2O')) = 40;
%         Vuptake(strcmpi(model.rxns,'exPI')) = 50;
    model.Vuptake = Vuptake;
end
    
switch lower(option)
    case 'fba'
        %calculate FBA fluxes and growth rate
        bounds.Vuptake = model.Vuptake;
        bounds.vl = zeros(model.nt_rxn,1);        
        bounds.vl(logical(model.rev)) = -100;
        bounds.vu = zeros(model.nt_rxn,1);          
        bounds.vu(bounds.vu==0) = 100;%bounds.Vuptake;
        [vLP,flag] = estimateLPgrowth(model,bounds);
        
        if flag
            %steady state FBA fluxes for direction
            model.Vss = vLP.v;
        else
            fprintf('FBA fluxes unassigned\n');
        end
    case 'pfba'
        %calcuate maximum objective (growth rate, etc)
        bounds.Vuptake = model.Vuptake;
        bounds.vl = zeros(model.nt_rxn,1);        
        bounds.vl(logical(model.rev)) = -100;
        bounds.vu = zeros(model.nt_rxn,1);          
        bounds.vu(bounds.vu==0) = 100;
        [vLPmax,~,model] = solveLP(model,bounds,model.bmrxn);
        model.vl(model.c==1) = -vLPmax.obj;
        
        %claculate pFBA fluxes and irreversible model
        model = run_pFBA(model);
end

    