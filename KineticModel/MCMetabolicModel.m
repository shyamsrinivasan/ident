function [allSolution,allfinalSS,FBAmodel,conc] =...
         MCMetabolicModel(FBAmodel,ensb,variable,initSolution)
if nargin < 4
    initSolution = struct([]);
end
nr(1) = FBAmodel.nt_metab;
nr(2) = FBAmodel.nint_metab;
nr(3) = FBAmodel.nt_rxn;
nr(4) = FBAmodel.n_rxn;
%Sample initial values and simulate ODEmodel
%Model Initialization
if isempty(initSolution)
    [FBAmodel,batch,solverP,saveData] = initializeModel(FBAmodel);
    %simulate models first to get initial SS
    %for isample = 1:nmodels
    %mname = sprintf('model%d',isample);
    mname = 'model1';
    [Solution,finalSS] =...
    callODEsolver(FBAmodel,ensb.(mname),variable,initSolution,batch,solverP);
%     %Obtaining final SS values for initial problem 
%     f_flag = 0;
%     i = size(Solution.y,2);
%     while ~f_flag
%     % for i = size(Solution.initSS.y,2):-1:1
%         if any(Solution.y(1:FBAmodel.nint_metab,i)<-solverP.MabsTol)
%             finalSS.y = Solution.y(:,i-1);
%             finalSS.t = Solution.t(i-1);
%             f_flag = 1;  
%         else
%             finalSS.y = Solution.y(:,i);
%             finalSS.t = Solution.t(i);
%             f_flag = 1;
%         end
%         i = i-1;
%     end
%     if ~f_flag
%         finalSS.t = Solution.t(end);
%         finalSS.y = Solution.y(:,end); 
%     end
    initSolution(1).(mname) = Solution;  
    %Calculate initial & final flux   
    initSolution.(mname).flux = calc_flux(FBAmodel,ensb.(mname),Solution.y(:,end));    
    allSolution.init = Solution;
    allfinalSS.init = finalSS;                
    %Calculate Flux
%     save_flux.t = intflux;                
%     allfinalSS.init.flux = intflux; 
    close all
end
%Monte Carlo Simulation of initial values
[FBAmodel,batch,solverP,saveData] = initializeModel(FBAmodel,10);
ninitsample = 5;%# samples
nvar = FBAmodel.nint_metab;
conc = zeros(nvar,ninitsample+1);
petconc = zeros(nvar,ninitsample);
flux = zeros(FBAmodel.nt_rxn,ninitsample+1);

%Perturb initial values
%Y0 range                
lb = [1e-4];
ub = [10];
pvar = {'C[c]'};
varname = {'A[c]','B[c]','C[c]','D[c]','E[c]','P[c]'};
% [y0new,pbind,savesimd] =...
% MCMetabolic_InitialVal(FBAmodel,pvar,finalSS.y,lb,ub,ninitsample);
[y0new,pbind,savesimd] =...
MCsample_initval(FBAmodel,[],pvar,finalSS.y,lb,ub,ninitsample);

for isamp = 1:ninitsample
    fprintf('\nStarting Sample %d of %d \n',isamp,ninitsample);                    
    ksampl = sprintf('sample_%d',isamp);
    [FBAmodel,batch,solverP,saveData] = initializeModel(FBAmodel);
    fprintf('Cell Growth Rate = %1.2g h-1\n',FBAmodel.gmax);
    saveData.filename = ksampl;
    %Choose a sampled value
    petconc(:,isamp) = y0new(:,isamp);
    %Save simulation conditions
    fname = sprintf('simconditions_%d',isamp);
    savefile(savesimd,fname,saveData);  
    %append pvar to  varname
%     for iv = 1:length(pvar)
%         if ~any(strcmpi(pvar{iv},varname))
%             try 
%                 varname = [varname;pvar{iv}];
%             catch
%                 varname = [varname pvar{iv}];
%             end
%         end
%     end
    %Add protein name
    pertbinitval.y = y0new(:,isamp);
    pertbinitval.pertbind = pbind;                              
    Solution = callODEsolver(FBAmodel,ensb.(mname),variable,pertbinitval,batch,solverP);    
    %Obtaining final SS values for initial problem 
    f_flag = 0;
    i = size(Solution.y,2);
    while ~f_flag
    % for i = size(Solution.initSS.y,2):-1:1
        if any(Solution.y(1:FBAmodel.nint_metab,i)<-solverP.MabsTol)
            finalSS.y = Solution.y(:,i-1);
            finalSS.t = Solution.t(i-1);
            f_flag = 1;  
        else
            finalSS.y = Solution.y(:,i);
            finalSS.t = Solution.t(i);
            f_flag = 1;
        end
        i = i-1;
    end
    if ~f_flag
        finalSS.t = Solution.t(end);
        finalSS.y = Solution.y(:,end); 
    end    
    allSolution.(ksampl) = Solution;
    allfinalSS.(ksampl) = finalSS;
%     savefile(Solution.initSS,ksampl,saveData);
    conc(:,isamp+1) = finalSS.y;
    %Calculate Flux
%     MC = finalSS.y(ng(1)+ng(2)+2:ng(1)+ng(2)+1+ng(3)+ng(4));
%     EC = finalSS.y(ng(1)+1:ng(1)+ng(2)+ng(3)+1)*1e-3;%umole to mmole
%     flux(:,isamp+1) = calc_flux(trnmodel,trnmodel.pmeter,MC,EC);
%     save_flux.y = flux(:,isamp+1);
%     allfinalSS.(ksampl).flux = flux(:,isamp+1);
    %Save Flux
%     fname = sprintf('flux_%d',isamp);
%     savefile(save_flux,fname,saveData);
    fprintf('Completed Sample #%d of %d\n',isamp,ninitsample);
    close all  
end 
conc(:,1) = allfinalSS.init.y;
% flux(:,1) = intflux;
close all
printMetResults(FBAmodel,allSolution,conc,petconc,[],varname);
%Plot Solution Curves
% [plotData,h_subfig] =...
% printResults(trnmodel,ng,allSolution,conc,petconc,flux,printvar,varname);
return

%Repeat MC Simulation
