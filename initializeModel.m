function [model,batch,solverP,saveData] = initializeModel(model,tmax,varargin)

if nargin < 3
    pmeter = [];
    variable = [];
else
    pmeter = varargin{1};
    variable = varargin{2};
end
if nargin < 2
    tmax = 500000;%s
end
batch = struct();
solverP = struct();

%Setting a few initial assignments and parameters that are unassigned
if ~isfield(model,'Vefflux')%Not needed ?
    model.Vefflux = 0;
end
%Kcat for all reactions
if ~isfield(model,'Kcat')
    model.Kcat = zeros(model.nt_rxn,1);
    model.Kcat(model.Kcat == 0) = 3000;%s-1
elseif isfield(model,'Kcat')
    if any(isnan(model.Kcat))
        model.Kcat(isnan(model.Kcat)) = 3000;
    end
    model.Kcat(model.Kcat == 0) = 3000;    
end
%Initial Extracellular concentrations
if ~isfield(batch,'init')
%     batch.init{1} = {'glc[e]';'lcts[e]';'ac[e]';'etoh[e]';'succ[e]';...
%                      'nh4[e]';'co2[e]';'acald[e]';'akg[e]';'for[e]';...
%                      'fum[e]';'h[e]';'h2o[e]';'o2[e]';'mal[e]'};
%     batch.init{2} = [2;0;0;0;0;0;0;0;0;0;0;0;0;1;0];%mmoles 
      batch.init{1} = {'A[e]','P[e]','D[e]','E[e]'};
      batch.init{2} = [20;0;0;0];
%       batch.init{1} = {'S[e]','B[e]','P[e]','A[e]'};
%       batch.init{2} = [20;0;0;0];
end
%??
if ~isfield(model,'kcat')
    model.kcat = 1;%s-1
end
if ~isfield(model,'ext_MC')
    model.ext_MC = assign_extconc(batch.init{1},batch.init{2},model);
end
%Fixed uptake rates for FBA
if ~isfield(model,'Vuptake')
    [~,Yflux] =...
    initializeConcentration(model,pmeter,variable,batch.init{1},...
                            batch.init{2},1);
%     nuprxns = length(model.Vupind);
%     model.Vuptake = zeros(nuprxns,1);
%     model.Vuptake(model.Vuptake==0) = 20;%mmole/gDCW.s
    model.Vuptake = 20;%Yflux(model.Vupind);
end

batch.tpmax = 10;%h  
%Integration Time
batch.tmax = tmax;%s
%ODE Solver parameters
solverP.RabsTol = 1e-6;
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
% ivec = 1;
% remVex = zeros(length(model.Vex),1);
% for irxn = 1:length(model.Vex)
%     nsubs = length(find(model.S(:,model.Vex(irxn))<0));
%     nprod = length(find(model.S(:,model.Vex(irxn))>0));
%     if nsubs > 1 && nprod > 1
%         model.V2rct(ivec) = model.Vex(irxn);
%         remVex(irxn) = 1;
%         ivec = ivec + 1;
%     end
% end
% if any(remVex)
%     model.Vex(logical(remVex)) = [];
% else
%     model.V2rct = [];
% end

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
[gMax,~,~,~,gMaxflag] = solveLP(model,'','',bounds,model.bmrxn);
% [pMax,~,~,~,pMaxflag] = solveLP(model,'','',bounds,find(strcmpi('Pex',model.rxns)));
fprintf('Uptake Flux = %2.3g\n',model.Vuptake);
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

return
%function calculateVariance()
%function sensAnalysis()