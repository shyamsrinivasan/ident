function model = estimateLPgrowth(model)

%Check if growth rate is possible
%Uptake Flux
bounds.Vuptake = model.Vuptake;
bounds.vl = zeros(model.nt_rxn,1);
% bounds.vl(bounds.vl==0) = -1;
bounds.vl(logical(model.rev)) = -100;%bounds.Vuptake;
bounds.vu = zeros(model.nt_rxn,1);          
%Corresponding flux bounds
bounds.vu(bounds.vu==0) = 100;%bounds.Vuptake;
%Determine Max and Min for flux to be constrained with =

[gMax,~,~,~,gMaxflag,~,model] = solveLP(model,'','',bounds,model.bmrxn);
% [pMax,~,~,~,pMaxflag] = solveLP(model,'','',bounds,find(strcmpi('Pex',model.rxns)));
%Print uptake fluxes
vupid = logical(model.Vuptake);
fprintf('Uptake Fluxes\n');
M = [model.rxns(vupid) num2cell(model.Vuptake(vupid))];
fprintf('%s \t %d\n',M{:});
% fprintf('Uptake Flux = %2.3g\n',model.Vuptake);
if ~isfield(model,'gmax') && gMaxflag > 0
    fprintf('Maximum feasible growth rate = %2.3g h-1\n',-gMax);
    model.gmax = 0.1;%-vMax;
elseif ~isfield(model,'gmax')
    model.gmax = 0.1;
elseif -gMax < model.gmax
    fprintf('Given maximum growth rate %2.3g is infeasible\n',model.gmax);
    fprintf('Maximum feasible growth rate = %2.3g h-1\n',-gMax);
    model.gmax = -gMax;
end  