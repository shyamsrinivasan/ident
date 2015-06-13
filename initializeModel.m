function [model,batch,solverP,saveData] = initializeModel(model,tmax)
if nargin < 2
    tmax = 500000;%s
end
batch = struct();
solverP = struct();
if ~isfield(model,'gmax')
    model.gmax = 0.05;%h-1
end
if ~isfield(model,'Vuptake')
    nuprxns = length(model.Vupind);
    model.Vuptake = zeros(nuprxns,1);
    model.Vuptake(model.Vuptake==0) = 1;%mmole/gDCW.s
end
if ~isfield(model,'Vefflux')
    model.Vefflux = 0;
end
if ~isfield(model,'Kcat')
    model.Kcat = zeros(model.nt_rxn,1);
    model.Kcat(model.Kcat == 0) = 3000;%s-1
elseif isfield(model,'Kcat')
    if any(isnan(model.Kcat))
        model.Kcat(isnan(model.Kcat)) = 3000;
    end
    model.Kcat(model.Kcat == 0) = 3000;    
end
if ~isfield(batch,'init')
    batch.init{1} = {'glc[e]';'lcts[e]';'ac[e]';'etoh[e]';'succ[e]';...
                     'nh4[e]';'co2[e]';'acald[e]';'akg[e]';'for[e]';...
                     'fum[e]';'h[e]';'h2o[e]';'o2[e]';'mal[e]';'pi[e]';'pyr[e]'};
    batch.init{2} = [200;0;0;0;0;0;0;0;0;0;0;2E-3;1000;1;0;10;0];%mmoles 
%       batch.init{1} = {'A[e]','D[e]','P[e]','E[e]'};
%       batch.init{2} = [1;0;0;0];
end
if ~isfield(model,'kcat')
    model.kcat = 1;%s-1
end
batch.tpmax = 10;%h  
%Integration Time
batch.tmax = tmax;%s
%ODE Solver parameters
solverP.RabsTol = 1e-7;
solverP.PabsTol = 1e-6;
solverP.MabsTol = 1e-5;
solverP.RelTol = 1e-4;
solverP.MaxIter = 1000;    
solverP.MaxDataPoints = 200;
solverP.tmax = tmax;%s
solverP.tout = 0.01;
%File save location/folder
saveData.filename = '';%sprintf('ExptCondition_%d',exptnum);
saveData.dirname =...
'C:\Users\shyam\Documents\Courses\CHE1125Project\Results\KModel';

%Assigning Indices
ivec = 1;
remVex = zeros(length(model.Vex),1);
for irxn = 1:length(model.Vex)
    nsubs = length(find(model.S(:,model.Vex(irxn))<0));
    nprod = length(find(model.S(:,model.Vex(irxn))>0));
    if nsubs > 1 && nprod > 1
        model.V2rct(ivec) = model.Vex(irxn);
        remVex(irxn) = 1;
        ivec = ivec + 1;
    end
end
if any(remVex)
    model.Vex(logical(remVex)) = [];
else
    model.V2rct = [];
end

%Check if growth rate is possible
%Uptake Flux
bounds.Vuptake = model.Vuptake;
bounds.vl = zeros(model.nt_rxn,1);
% bounds.vl(bounds.vl==0) = -1;
bounds.vl(logical(model.reversible)) = -100;%bounds.Vuptake;
bounds.vu = zeros(model.nt_rxn,1);          
%Corresponding flux bounds
bounds.vu(bounds.vu==0) = 100;%bounds.Vuptake;
%Determine Max and Min for flux to be constrained with =
[vMax,~,~,~,Maxflag] = solveLP(model,'','',bounds,model.bmrxn);
if Maxflag > 0 && -vMax < model.gmax
    fprintf('Given Maximum growth rate %2.3g is infeasible\n',model.gmax);
    model.gmax = -vMax;
end
return
%function calculateVariance()
%function sensAnalysis()