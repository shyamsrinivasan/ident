function [allSolution,allfinalSS] =...
         MCsimulation(model,parameter,variable,nsamples,pvar,lb,ub,...
                      allSolution,allfinalSS,varname)

nvar = model.nt_metab;
conc = zeros(nvar,nsamples+1);
petconc = zeros(nvar,nsamples);
flux = zeros(model.nt_rxn,nsamples+1);

[y0new,pbind,savesimd] =...
MCsample_initval(model,[],pvar,allfinalSS.init.y,lb,ub,nsamples);

for isamp = 1:nsamples
    fprintf('\nStarting Sample %d of %d \n',isamp,nsamples);                    
    ksp_name = sprintf('sample_%d',isamp);
    [model,batch,solverP,saveData] = initializeModel(model,300);
    fprintf('Cell Growth Rate = %1.2g h-1\n',model.gmax);
    saveData.filename = ksp_name;
    
    %Choose a sampled value
    petconc(:,isamp) = y0new(:,isamp);
    
    %Save simulation conditions
    fname = sprintf('simconditions_%d',isamp);
    savefile(savesimd,fname,saveData);  

    %Add protein name
    pertbinitval.y = y0new(:,isamp);
    pertbinitval.pertbind = pbind;                              
    [Solution,finalSS] =...
    callODEsolver(model,parameter,variable,pertbinitval,batch,solverP);    
    
    allSolution.(ksp_name) = Solution;
    allfinalSS.(ksp_name) = finalSS;
    savefile(Solution,ksp_name,saveData);
    conc(:,isamp+1) = finalSS.y;
    
    %Calculate Flux
    flux(:,isamp+1) = calc_flux(model,parameter,finalSS.y);
    save_flux.t = [];
    save_flux.y = flux(:,isamp+1);
    allfinalSS.(ksp_name).flux = flux(:,isamp+1);
    
    %Save Flux
    fname = sprintf('flux_%d',isamp);
    savefile(save_flux,fname,saveData);
    fprintf('Completed Sample #%d of %d\n',isamp,nsamples);
    close all  
end 
conc(:,1) = allfinalSS.init.y;
flux(:,1) = allfinalSS.init.flux;
close all

%Plot Solution Curves
[hfig,hsubfig] =...
printMetResults(model,allSolution,conc,petconc,[],varname);

%Plot Fluxes & Bin and plot Flux Distribution
printvar = {'Pin','P1','P2','P3','P4','P5','P6','P7','P8','P9','Biomass'};
plotflux_bar(model,flux./model.Vuptake,printvar);

%Plot Steady State Concentrations
plotSSexpression(model,[],conc,petconc,varname);

%Calculate & Plot Envelope
[hsubfig2,prxnid,flag] = FluxEnvelope(model);

%Plot Scatter within Envelope (or superimpose)
if flag > 0
    plotflux_envelope(model,flux./model.Vuptake,printvar,hsubfig2,prxnid);
else
    fprintf('\n Envelopes Infeasible for growth rate = %3.2g h-1\n',model.gmax);
end  

%Plot Flux Scatter
plotflux_scatter(model,flux./model.Vuptake,printvar);

%Plot Solution Curves
% [plotData,h_subfig] =...
% printResults(trnmodel,ng,allSolution,conc,petconc,flux,printvar,varname);
return