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
if nargin<2
    option = 'any';
end

%fix flux uptakes for FBA solution
Vuptake = fixUptake(model,Vup_struct);
if ~isfield(model,'Vuptake')
    model.Vuptake = Vuptake;
end
% bounds = changebounds(model);

% change reaction flux bounds to desired values
[model,bounds] = changebounds(model,ess_rxn);
    
switch lower(option)
    case 'fba'
        %calculate FBA fluxes and growth rate
        model = calculateFBAfluxes(model);
        
    case 'pfba'       
        model = calculatePFBAfluxes(model);
        
    case 'any'
        % calculate fluxes using either fba or pfba whichever is feasible        
        % first use pfba
        try
            model = calculatePFBAfluxes(model);
        catch % calculate pfba
            fprintf('\n pFBA was unsuccesful \n Trying FBA now...\n');
            model = calculateFBAfluxes(model);
        end        
    otherwise
        fprintf('Not calculating FBA SS fluxes. \nOnly Calculating bounds.\n');
        [~,~,model] = solveLP(model,bounds,ess_rxn,model.bmrxn,Vup_struct);
end

function model = calculatePFBAfluxes(model)
% calcuate maximum objective (growth rate, etc)
if ~prxnid
    prxnid = model.bmrxn;
end
% atp = find(strcmpi(model.rxns,'ATPM'));
[vLPmax,~,model] = solveLP(model,bounds,ess_rxn,prxnid,Vup_struct);
model.vl(model.c==1) = -vLPmax.obj;

% claculate pFBA fluxes and irreversible model
vLPmin = run_pFBA(model,ess_rxn,Vup_struct);

if vLPmin.flag
    model.Vss = vLPmin.v;
end
fprintf('pFBA was successful\n');
end

function model = calculateFBAfluxes(model)
% Determine Max and Min for flux to be constrained with
if ~prxnid
    prxnid = model.bmrxn;
end        
[vLPmax,vLPmin,model] = solveLP(model,bounds,ess_rxn,prxnid,Vup_struct);

if vLPmax.flag
    % steady state FBA fluxes for direction
    model.Vss = vLPmax.v;
    fprintf('FBA successful\n');
else
    fprintf('FBA fluxes unassigned\n');
end
end

end  