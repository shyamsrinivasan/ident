function model = FBAfluxes(model,option,ess_rxn,Vup_struct)
if nargin <4
    Vup_struct = ([]);
end

%fix flux uptakes for FBA solution
model = fixUptake(model,Vup_struct);
    
switch lower(option)
    case 'fba'
        %calculate FBA fluxes and growth rate
        bounds.Vuptake = model.Vuptake;
        bounds.vl = zeros(model.nt_rxn,1);        
        bounds.vl(logical(model.rev)) = -100;
        bounds.vu = zeros(model.nt_rxn,1);          
        bounds.vu(bounds.vu==0) = 100;%bounds.Vuptake;
%         [vLP,flag] = estimateLPgrowth(model,ess_rxn,bounds);
        %Determine Max and Min for flux to be constrained with
        [vLPmax,vLPmin,model] = solveLP(model,ess_rxn,bounds,model.bmrxn,Vup_struct);
        
        if flag
            %steady state FBA fluxes for direction
            model.Vss = vLPmax.v;
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
%         atp = find(strcmpi(model.rxns,'ATPM'));
        [vLPmax,~,model] = solveLP(model,ess_rxn,bounds,model.bmrxn,Vup_struct);
        model.vl(model.c==1) = -vLPmax.obj;
        
        %claculate pFBA fluxes and irreversible model
        model = run_pFBA(model,ess_rxn,Vup_struct);
    otherwise
        fprintf('Not calculating FBA SS fluxes. \nOnly Calculating bounds.\n');
        bounds.Vuptake = model.Vuptake;
        bounds.vl = zeros(model.nt_rxn,1);        
        bounds.vl(logical(model.rev)) = -100;
        bounds.vu = zeros(model.nt_rxn,1);          
        bounds.vu(bounds.vu==0) = 100;%bounds.Vuptake;
        [~,~,model] = solveLP(model,ess_rxn,bounds,model.bmrxn,Vup_struct);
end

    