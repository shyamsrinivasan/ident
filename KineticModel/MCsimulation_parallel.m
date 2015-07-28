function [allSolution,allfinalSS,varargout] =...
         MCsimulation_parallel(model,parameter,variable,nsamples,pvar,lb,ub,...
                      allSolution,allfinalSS,varname)

%Parallel Implementation
% cluster = parcluster;
if isempty(gcp('nocreate'))
    parpool(4);
end

nvar = model.nt_metab;
conc = zeros(nvar,nsamples+1);
petconc = zeros(nvar,nsamples);
flux = zeros(model.nt_rxn,nsamples+1);
petflux = zeros(model.nt_rxn,nsamples);
if ~isempty(pvar)
    if any(strcmpi(pvar{1},model.mets))
        [y0new,pbind,savesimd,status] =...
        MCsample_initval(model,[],pvar,allfinalSS.init.y,lb,ub,nsamples);
    else
        %Do Vmax perturbation simulation
        %Get samples for Vmax
    end
end
if status >= 0
    samp_name = cell(nsamples,1);
    for isamp = 1:nsamples
        ksp_name = sprintf('sample_%d',isamp);
        samp_name{isamp} = ksp_name;
        [model,batch,solverP,saveData] = initializeModel(model,300000);
        data.(ksp_name).model = model;
        data.(ksp_name).batch = batch;
        data.(ksp_name).solverP = solverP;
        data.(ksp_name).saveData = saveData;
    end
    ParSol = cell(nsamples,1);
    ParfSS = cell(nsamples,1);
else
    return
end

parfor isamp = 1:nsamples
    fprintf('\nStarting Sample %d of %d \n',isamp,nsamples);                    
%     ksp_name = sprintf('sample_%d',isamp);
    saveData = data.(samp_name{isamp}).saveData;
    fprintf('Cell Growth Rate = %1.2g h-1\n',data.(samp_name{isamp}).model.gmax);
    saveData.filename = samp_name{isamp};
    
    %Choose a sampled value
    petconc(:,isamp) = y0new(:,isamp);
    petflux(:,isamp) = calc_flux(model,parameter,y0new(:,isamp));
    
    %Save simulation conditions
    fname = sprintf('simconditions_%d',isamp);
%     savefile(savesimd,fname,saveData);  
    
    %Create job
%     job.(ksp_name) = createJob(cluster);
    
    %Add protein name
    pertbinitval = struct();
    pertbinitval.y = y0new(:,isamp);
    pertbinitval.flux = petflux(:,isamp);
    pertbinitval.pertbind = pbind;  
%     createTask(job.(ksp_name),@callODEsolver,3,...
%     {model,parameter,variable,pertbinitval,batch,solverP});    
    [Solution,finalSS] =...
    callODEsolver(data.(samp_name{isamp}).model,...
    parameter,variable,pertbinitval,data.(samp_name{isamp}).batch,...
    data.(samp_name{isamp}).solverP);  
    
    ParSol{isamp} = Solution;
    ParfSS{isamp} = finalSS;
%     savefile(ParSol{isamp},samp_name{isamp},saveData);
    conc(:,isamp+1) = ParfSS{isamp}.y;
    
    %Calculate Flux
    flux(:,isamp+1) = calc_flux(model,parameter,ParfSS{isamp}.y);   

end
% delete(gcp);
% for isamp = 1:nsamples
    
% submit(j);
% wait(j);
% ParallelData = fetchOutputs(j);
% delete(j);
fprintf('Saving data to file ... \n');
for isamp = 1:nsamples
%     fprintf('\nStarting Sample %d of %d \n',isamp,nsamples);                    
%     ksp_name = sprintf('sample_%d',isamp);
    
    model = data.(samp_name{isamp}).model;
    saveData = data.(samp_name{isamp}).saveData;
    saveData.filename = samp_name{isamp};
    
    allSolution.(samp_name{isamp}) = ParSol{isamp};%ParallelData{isamp,1};
    allfinalSS.(samp_name{isamp}) = ParfSS{isamp};%ParallelData{isamp,2};
%     savefile(allSolution.(samp_name{isamp}),samp_name{isamp},saveData);
%     conc(:,isamp+1) = allfinalSS.(ksp_name).y;
    
    %Calculate Flux
%     flux(:,isamp+1) = calc_flux(model,parameter,allfinalSS.(ksp_name).y);
    save_flux.t = [];
    save_flux.y = flux(:,isamp+1);
    allfinalSS.(samp_name{isamp}).flux = flux(:,isamp+1);
    
    %Save Flux
    fname = sprintf('flux_%d',isamp);
%     savefile(save_flux,fname,saveData);
%     fprintf('Completed Sample #%d of %d\n',isamp,nsamples);
    close all  
end
fprintf('Complete\n');
conc(:,1) = allfinalSS.init.y;
flux(:,1) = allfinalSS.init.flux;
close all

%Plot Solution Curves
% [hfig,hsubfig] =...
% printMetResults(model,allSolution,conc,petconc,[],varname);

%Bin Concentrations/Fluxes
[conc,MSSconc] = binConcentrations(conc);
[flux,MSSflux] = binConcentrations(flux);

[petconc,MSSpconc] = binConcentrations(petconc);
[petflux,MSSpflux] = binConcentrations(petflux);

%Plot Fluxes & Bin and plot Flux Distribution
% printvar = {'Pin','v1','v2','v3','v4','v5','v6','v7','v8','Pout','Bout','Aout','BiomassEX'};
printvar = {'Pin','v1','v2','v3','v4','v5','v6','Pout','Dout','Eout','BiomassEX'};
plotflux_bar(model,flux,printvar);

%Plot Steady State Concentrations
plotSSexpression(model,[],conc,petconc,varname,'concentration');

%Plot Steady State Fluxes
plotSSexpression(model,[],flux,petflux,printvar,'flux');

%Calculate & Plot Envelope
prdid = strcmpi('Pex',model.rxns);
[hsubfig2,prxnid,flag] = FluxEnvelope(model,printvar,prdid);

%Plot Scatter within Envelope (or superimpose)
% notbm = setdiff(1:size(flux,1),model.bmrxn);
if flag > 0
    plotflux_envelope(model,flux,printvar,hsubfig2,prxnid);
%     plotflux_envelope(model,petflux,printvar,hsubfig2,prxnid,2);
else
    fprintf('\n Envelopes Infeasible for growth rate = %3.2g h-1\n',model.gmax);
end   
% getFigure(petconc',petflux',model);


Csigma = OutputVariability(conc);
Vsigma = OutputVariability(flux);

CinSigma = OutputVariability(petconc);
Vinsigma = OutputVariability(petflux);

varargout{1} = y0new;
varargout{2} = MSSconc;
varargout{3} = MSSpconc;
varargout{4} = MSSflux;
varargout{5} = MSSpflux;


%Plot Flux Scatter
% plotflux_scatter(model,flux(notbm,:)./model.Vuptake,printvar);

%Plot Distribution of fluxes & concentrations (initial and final)
% plotfluxdistr(model,flux,printvar)


%Plot Solution Curves
% [plotData,h_subfig] =...
% printResults(model,allSolution,[],[],[],varname);

% printMetResults(model,allSolution,[],[],[],varname);

return
