function [model,solverP,saveData] = imodel(model,type,tmax,ess_rxn,Vup_struct)
if nargin <5
    Vup_struct = ([]);
end
if nargin<4
    ess_rxn = {};
end
if nargin<3
    tmax = 5000;%s
end
solverP = struct();

switch type
    case 'ode'
        % ODE Solver parameters
        solverP.RabsTol = 1e-12;
        solverP.PabsTol = 1e-12;
        solverP.MabsTol = 1e-12;
        solverP.RelTol = 1e-14;
        solverP.MaxIter = 70000;    
        solverP.MaxDataPoints = 500;
        solverP.tmax = tmax;%s
        solverP.tout = 0.01;
    case 'nla'
        solverP.FuncNormTol = 1e-10;
        solverP.ScaledStepTol = 1e-5;
        solverP.KrylovMaxDim = 10;
        solverP.MaxNumRestarts = 2;
        solverP.MaxNumSetups = 5;
        solverP.MaxNumSteps = 500;
        solverP.Verbose = true;
end

% data file save location/folder
saveData.filename = '';%sprintf('ExptCondition_%d',exptnum);
saveData.dirname =...
'C:\Users\shyam\Documents\Courses\CHE1125Project\Results\KModel';

%--------------------------------------------------------------------------
% recheck growth rates
% fix flux uptakes for FBA solution
Vuptake = fixUptake(model,Vup_struct);
if ~isfield(model,'Vuptake')
    model.Vuptake = Vuptake;
end
% change reaction flux bounds to desired values
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
