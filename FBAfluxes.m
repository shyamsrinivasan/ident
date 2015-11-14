function model = FBAfluxes(model,option,ess_rxn,Vup_struct,prxnid)
if nargin<5 || isempty(prxnid)
    prxnid = 0;
end
if nargin <4
    Vup_struct = ([]);
end
if nargin< 3
    ess_rxn = {};
end

%fix flux uptakes for FBA solution
model = fixUptake(model,Vup_struct);
% bounds = changebounds(model);
% 
% bounds.Vuptake = model.Vuptake;
% bounds.vl = zeros(model.nt_rxn,1);        
% bounds.vl(logical(model.rev)) = -100;
% bounds.vu = zeros(model.nt_rxn,1);          
% bounds.vu(bounds.vu==0) = 100;%bounds.Vuptake;

[model,bounds] = changebounds(model,ess_rxn);
    
switch lower(option)
    case 'fba'
        %calculate FBA fluxes and growth rate
%         bounds.Vuptake = model.Vuptake;
%         bounds.vl = zeros(model.nt_rxn,1);        
%         bounds.vl(logical(model.rev)) = -100;
%         bounds.vu = zeros(model.nt_rxn,1);          
%         bounds.vu(bounds.vu==0) = 100;%bounds.Vuptake;
%         [vLP,flag] = estimateLPgrowth(model,ess_rxn,bounds);
        %Determine Max and Min for flux to be constrained with
        if ~prxnid
            prxnid = model.bmrxn;
        end        
        [vLPmax,vLPmin,model] = solveLP(model,bounds,ess_rxn,prxnid,Vup_struct);
        
        if vLPmax.flag
            %steady state FBA fluxes for direction
            model.Vss = vLPmax.v;
        else
            fprintf('FBA fluxes unassigned\n');
        end
    case 'pfba'
        %calcuate maximum objective (growth rate, etc)
%         bounds.Vuptake = model.Vuptake;
%         bounds.vl = zeros(model.nt_rxn,1);        
%         bounds.vl(logical(model.rev)) = -100;
%         bounds.vu = zeros(model.nt_rxn,1);          
%         bounds.vu(bounds.vu==0) = 100;
        if ~prxnid
            prxnid = model.bmrxn;
        end
%         atp = find(strcmpi(model.rxns,'ATPM'));
        [vLPmax,~,model] = solveLP(model,bounds,ess_rxn,prxnid,Vup_struct);
        model.vl(model.c==1) = -vLPmax.obj;
        
        %claculate pFBA fluxes and irreversible model
        model = run_pFBA(model,ess_rxn,Vup_struct);
    otherwise
        fprintf('Not calculating FBA SS fluxes. \nOnly Calculating bounds.\n');
%         bounds.Vuptake = model.Vuptake;
%         bounds.vl = zeros(model.nt_rxn,1);        
%         bounds.vl(logical(model.rev)) = -100;
%         bounds.vu = zeros(model.nt_rxn,1);          
%         bounds.vu(bounds.vu==0) = 100;%bounds.Vuptake;
        [~,~,model] = solveLP(model,bounds,ess_rxn,model.bmrxn,Vup_struct);
end

    