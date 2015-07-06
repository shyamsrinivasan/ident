function [inSolution,infinalSS,varargout] =...
MCmodel_parallel(model,inSolution,ensb,variable,nmodels,varname)

%Parallel Implementation
if isempty(gcp('nocreate'))
    parpool(4);
end

nvar = model.nt_metab;
conc = zeros(nvar,nmodels);
flux = zeros(model.nt_rxn,nmodels);
model_name = cell(nmodels,1);

%Initializzation Loop
for imodel = 1:nmodels
    mname = sprintf('model%d',imodel);    
    fprintf('%s\n',mname);
    model_name{imodel} = mname;
    [model,batch,solverP,saveData] = initializeModel(model,50,ensb.(mname),variable);
    %model
    data.(mname).model = model;
    %Inputs
    data.(mname).batch = batch;
    %solver Options
    data.(mname).solverP = solverP;
    data.(mname).saveData = saveData;
    %Empty initial values (assigned before simulation in callODEsolver
    if ~isfield(inSolution,mname)
%     if isempty(inSolution)
        inSolution.(mname) = struct([]);
    end   
end
ModSS = cell(nmodels,1);
ModFSS = cell(nmodels,1);

%parallel loop using parfor
parfor imodel = 1:nmodels
    fprintf('\nStarting Model %d of %d \n',imodel,nmodels);   
    saveData = data.(model_name{imodel}).saveData;
    fprintf('Cell Growth Rate = %1.2g h-1\n',data.(model_name{imodel}).model.gmax);
    saveData.filename = model_name{imodel};
    
    %Save simulation conditions
    fname = sprintf('model_%d',imodel);
%     savefile(savesimd,fname,saveData); 
    
    %simulate models to get SS for given model
    [Solution,finalSS] =...
    callODEsolver(data.(model_name{imodel}).model,...
                  ensb.(mname),variable,inSolution.(mname),...
                  data.(model_name{imodel}).batch,...
                  data.(model_name{imodel}).solverP);
    ModSS{imodel} = Solution;
    ModFSS{imodel} = finalSS;
%     savefile(ModSS{imodel},model_name{imodel},saveData);
    conc(:,imodel) = finalSS.y;
    
    %Calculate Flux
    flux(:,imodel) = calc_flux(model,ensb.(mname),finalSS.y);  
    ModFSS{imodel}.flux = flux(:,imodel);
end

fprintf('Saving data to file ... \n');
for imodel = 1:nmodels
%     fprintf('\nStarting Sample %d of %d \n',isamp,nsamples);                    
%     ksp_name = sprintf('sample_%d',isamp);
    
    model = data.(model_name{imodel}).model;
    saveData = data.(model_name{imodel}).saveData;
    saveData.filename = model_name{imodel};
    
    inSolution.(model_name{imodel}) = ModSS{imodel};%ParallelData{isamp,1};
    infinalSS.(model_name{imodel}) = ModFSS{imodel};%ParallelData{isamp,2};
%     savefile(allSolution.(samp_name{isamp}),samp_name{isamp},saveData);
%     conc(:,isamp+1) = allfinalSS.(ksp_name).y;
    
    
    %Calculate Flux
    save_flux.t = [];
    save_flux.y = flux(:,imodel);
%     infinalSS.(model_name{imodel}).flux = flux(:,imodel);
    
    %Save Flux
%     fname = sprintf('flux_%d',isamp);
%     savefile(save_flux,fname,saveData);
%     fprintf('Completed Sample #%d of %d\n',isamp,nsamples);     
end
fprintf('Complete\n');
close all 

%Plot Fluxes & Bin and plot Flux Distribution
printvar = {'Pin','P1','P2','P3','P4','P5','P6','P7','P8','P9','BiomassEX'};
plotflux_bar(model,flux,printvar);

%Plot Steady State Concentrations
plotSSexpression(model,[],conc,[],varname,'concentration');

%Calculate & Plot Envelope
[hsubfig2,prxnid,flag] = FluxEnvelope(model,printvar);

%Plot Scatter within Envelope (or superimpose)
% notbm = setdiff(1:size(flux,1),model.bmrxn);
if flag > 0
    plotflux_envelope(model,flux,printvar,hsubfig2,prxnid);
%     plotflux_envelope(model,petflux,printvar,hsubfig2,prxnid,2);
else
    fprintf('\n Envelopes Infeasible for growth rate = %3.2g h-1\n',model.gmax);
end  

varargout{1} = conc;
varargout{2} = flux;






