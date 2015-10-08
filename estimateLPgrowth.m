function [vLPmax,flag,vLPmin,model] = estimateLPgrowth(model,bounds)
if nargin <2
    %Check if growth rate is possible
    %Uptake Flux
    bounds.Vuptake = model.Vuptake;
    bounds.vl = zeros(model.nt_rxn,1);
    % bounds.vl(bounds.vl==0) = -1;
    bounds.vl(logical(model.rev)) = -100;%bounds.Vuptake;
    bounds.vu = zeros(model.nt_rxn,1);          
    %Corresponding flux bounds
    bounds.vu(bounds.vu==0) = 100;%bounds.Vuptake;
end

flag = 1;

%Determine Max and Min for flux to be constrained with
[vLPmax,vLPmin,model] = solveLP(model,bounds,model.bmrxn);

%Print uptake fluxes
vupid = logical(model.Vuptake);
fprintf('Uptake Fluxes\n');
M = [model.rxns(vupid) num2cell(model.Vuptake(vupid))];
fprintf('%s \t %d\n',M{:});

%print maximization result
if vLPmax.flag > 0
    fprintf('Maximum feasible product = %2.3g \n',-vLPmax.obj);
else
    fprintf('Maximization Infeasible\n');
    flag = -1;
%     model.gmax = 0.1;%-vMax;
% elseif ~isfield(model,'gmax')
%     model.gmax = 0.1;
% elseif -gMax < model.gmax
%     fprintf('Given maximum growth rate %2.3g is infeasible\n',model.gmax);
%     fprintf('Maximum feasible growth rate = %2.3g h-1\n',-gMax);
%     model.gmax = -gMax;
end  

%print minization result
if vLPmin.flag > 0
    fprintf('Minimum feasible product = %2.3g h-1\n',vLPmin.obj);
else
    fprintf('Minimization Infeasible\n');
    if flag < 0
        flag = -3;
    else
        flag = -2;
    end
end