function [vLPmax,flag,vLPmin,model] = estimateLPgrowth(model,ess_rxn,Vup_struct,bounds)

if nargin <4
    %Check if growth rate is possible
    %Uptake Flux
%     bounds.Vuptake = model.Vuptake;
%     bounds.vl = zeros(model.nt_rxn,1);
%     bounds.vl(logical(model.rev)) = -100;%bounds.Vuptake;
%     bounds.vu = zeros(model.nt_rxn,1);          
%     bounds.vu(bounds.vu==0) = 100;%bounds.Vuptake;
    [model,bounds] = changebounds(model,ess_rxn);
end
if nargin<3
    Vup_struct([]);
end
if nargin<2
    ess_rxn = {};
end

flag = 1;

%Determine Max and Min for flux to be constrained with using FBA
[vLPmax,vLPmin] = solveLP(model,bounds,ess_rxn,model.bmrxn);

%Determine Max and Min for flux to be constrained with using pFBA
% model.vl(model.c==1) = -vLPmax.obj;
% [~,Vss] = run_pFBA(model,ess_rxn,Vup_struct);

if vLPmax.flag > 0
    if abs((-vLPmax.obj)-model.Vss(model.bmrxn))<=1e-8
        maxGr = -vLPmax.obj;
    else
        warning('estLP:muDiff',...
        'FBA and pFBA growth rates are different\nThe process will proceed regardless with pFBA growth rates');
        maxGr = model.Vss(model.bmrxn);
    end
else
    error('estLP:LPfeas','The FBA problem is infeasible');
end

%Print uptake fluxes
vupid = logical(model.Vuptake);
fprintf('Uptake Fluxes\n');
M = [model.rxns(vupid) num2cell(model.Vuptake(vupid))];
fprintf('%s \t %d\n',M{:});

%print maximization result
if vLPmax.flag > 0
    fprintf('Maximum feasible product = %2.3g \n',maxGr);
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