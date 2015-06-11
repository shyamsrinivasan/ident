function [allSolution,conc,petconc,flux] =...
         MCdynamic(trnmodel,FBAmodel,defparval,ng,varname,ndynsample,initSolution)
if isempty(initSolution)
    [trnmodel,batch,solverP,saveData] = initializeModel(trnmodel);
    %calculate solution for model at gmax and set initSolution 
    [~,Solution,~,finalSS] =...
    trnperturbation(trnmodel,FBAmodel,batch,defparval,ng,...
                    solverP,initSolution,...
                    varname,saveData); 
    initSolution = Solution.initSS;
    allSolution.init = Solution.initSS;
    allfinalSS.init = finalSS;                
    %Calculate Flux
    MC = finalSS.y(ng(1)+ng(2)+2:ng(1)+ng(2)+1+ng(3)+ng(4));
    EC = finalSS.y(ng(1)+1:ng(1)+ng(2)+ng(3)+1)*1e-3;%umole to mmole
    intflux = calc_flux(trnmodel,trnmodel.pmeter,MC,EC);
    save_flux.t = intflux;                
    allfinalSS.init.flux = intflux;                
end

conc = zeros(sum(ng)+1,ndynsample+1);
petconc = zeros(sum(ng)+1,ndynsample);
flux = zeros(trnmodel.nt_rxn,ndynsample+1);
for isamp = 1:ndynsample
    ksampl = sprintf('sample_%d',isamp);
    [trnmodel,batch,solverP,saveData] = initializeModel(trnmodel);
    saveData.filename = ksampl;
    %Perturb initial values
    y0 = allfinalSS.init.y;
    %Y0 range
    lb = [1 1];
    ub = [100 100];
    pvar = {'P6','P5'};
    [y0new,pbind,savesimd] = MCsample_initval(trnmodel,ng,pvar,y0,lb,ub);
    petconc(:,isamp) = y0new;
    %Save simulation conditions
    fname = sprintf('simconditions_%d',isamp);
    savefile(savesimd,fname,saveData);  
    %append pvar to  varname
    for iv = 1:length(pvar)
        if ~any(strcmpi(pvar{iv},varname))
            try 
                varname = [varname;pvar{iv}];
            catch
                varname = [varname pvar{iv}];
            end
        end
    end
    %Add protein name
    pertbinitval.y = y0new;
    pertbinitval.pertbind = pbind;                       
    trnmodel.ext_MC = [];
    [~,Solution,~,finalSS] =...
    trnperturbation(trnmodel,FBAmodel,batch,defparval,ng,...
                    solverP,pertbinitval,...
                    varname,saveData); 
%                 initSolution = Solution.initSS;
    allSolution.(ksampl) = Solution.initSS;
    allfinalSS.(ksampl) = finalSS;
    savefile(Solution.initSS,ksampl,saveData);
    conc(:,isamp+1) = finalSS.y;
    %Calculate Flux
    MC = finalSS.y(ng(1)+ng(2)+2:ng(1)+ng(2)+1+ng(3)+ng(4));
    EC = finalSS.y(ng(1)+1:ng(1)+ng(2)+ng(3)+1)*1e-3;%umole to mmole
    flux(:,isamp+1) = calc_flux(trnmodel,trnmodel.pmeter,MC,EC);
    save_flux.y = flux(:,isamp+1);
    allfinalSS.(ksampl).flux = flux(:,isamp+1);
    %Save Flux
    fname = sprintf('flux_%d',isamp);
    savefile(save_flux,fname,saveData);
    fprintf('Completed Sample #%d of %d\n',isamp,ndynsample);
    close all  
end  
conc(:,1) = allfinalSS.init.y;
flux(:,1) = intflux;
return