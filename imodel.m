function [model,solverP,saveData] = imodel(model,tmax,ess_rxn,Vup_struct)
if nargin <4
    Vup_struct = ([]);
end
if nargin< 3
    ess_rxn = {};
end
if nargin < 2
    tmax = 500000;%s
end
solverP = struct();

%ODE Solver parameters
solverP.RabsTol = 1e-10;
solverP.PabsTol = 1e-10;
solverP.MabsTol = 1e-10;
solverP.RelTol = 1e-12;
solverP.MaxIter = 2000;    
solverP.MaxDataPoints = 1000;
solverP.tmax = tmax;%s
solverP.tout = 0.01;

%data file save location/folder
saveData.filename = '';%sprintf('ExptCondition_%d',exptnum);
saveData.dirname =...
'C:\Users\shyam\Documents\Courses\CHE1125Project\Results\KModel';

%--------------------------------------------------------------------------
%recheck growth rates
%fix flux uptakes for FBA solution
Vuptake = fixUptake(model,Vup_struct);
if ~isfield(model,'Vuptake')
    model.Vuptake = Vuptake;
end
%change reaction flux bounds to desired values
[model,bounds] = changebounds(model,ess_rxn);
if isfield(model,'bmrxn')
    prxnid = model.bmrxn;
end  

if ~isempty(prxnid)
    vLPmax = solveLP(model,bounds,ess_rxn,prxnid,Vup_struct);
    % [vLPmax,~,~,model] = estimateLPgrowth(model,ess_rxn,Vup_struct);
    if vLPmax.flag
        fprintf('Maximum feasible growth rate = %2.3g h-1\n',-vLPmax.obj);
    else
        fprintf('Growth is Infeasible\n');
    end
else
    fprintf('Model does not have a growth equation.\nContinuing...\n');
end
