function [model,solverP,saveData] = imodel(model,ess_rxn,Vup_struct,tmax)
if nargin < 4
    tmax = 500000;%s
end
solverP = struct();

%ODE Solver parameters
solverP.RabsTol = 1e-8;
solverP.PabsTol = 1e-8;
solverP.MabsTol = 1e-8;
solverP.RelTol = 1e-10;
solverP.MaxIter = 1000;    
solverP.MaxDataPoints = 200;
solverP.tmax = tmax;%s
solverP.tout = 0.01;

%data file save location/folder
saveData.filename = '';%sprintf('ExptCondition_%d',exptnum);
saveData.dirname =...
'C:\Users\shyam\Documents\Courses\CHE1125Project\Results\KModel';

%recheck growth rate
[vLPmax,~,~,model] = estimateLPgrowth(model,ess_rxn,Vup_struct);
if vLPmax.flag
    fprintf('Maximum feasible growth rate = %2.3g h-1\n',-vLPmax.obj);
else
    fprintf('Growth is Infeasible\n');
end
