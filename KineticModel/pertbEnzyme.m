% function [Sol] = pertbEnzyme(enz_struct,model,pmeter,variable,inSol,SolverOpt)
%enzyme perturbations
%Simulate every model in the ensemble for every Vmax sample
%Inside 'for loop of samples'
function [Sol] = pertbEnzyme(enz_struct,model,pmeter,variable,inSol,SolverOpt)

enz_tf = strcmpi(enz_struct.enzname,model.Enzyme);
if any(enz_tf)
    Vmax = pmeter.Vmax(enz_tf);
    dVmax = Vmax*enz_struct.percent/100;%percent% change in Vmax
    %Increase Vmax
    if strcmpi(enz_struct.change,'increase')
        pmeter.Vmax(enz_tf) = Vmax + dVmax;
        Sol = callODEsolver(model,pmeter,variable,inSol,SolverOpt);
    %Decrease Vmax
    elseif strcmpi(enz_struct.change,'decrease')
        pmeter.Vmax(enz_tf) = Vmax - dVmax;
        Sol = callODEsolver(model,pmeter,variable,inSol,SolverOpt);
    end
else
    fprintf('No Matching Enzyme\n');
end
% close all
return