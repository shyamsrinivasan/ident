% function [initval,Solution,status,finalSS,flag,indx,...
%           hallfig,defparval,exptnum] = trnperturbation(model,FBAmodel,batch,...
%           defparval,ng,solverP,initSolution,varname,saveData,exptnum)
%**************************************************************************
% Function to calculate initial SS solution and the effect of perturbations 
% to inputs
% initval
% Solution - Solution structure with time(t) and variable(y) fields
% status
% finalSS - Final SS value taken as the Solution.initSS.y(:end) 
% flag
% indx
% hallfig - Figure objects for all figures drawn within the function by 
%           calling dynamicplot
% defparval - Default parameter values
% exptnum
% ******Inputs
% model - Strructure of the integrated model
% FBAmodel - Structure of the metabolic model
% batch - Structure of input information to simulate model
%         Fields also include maximum simulation time tmax
% defparval - Structure of default parameter values. Structure fields
%             correspond to parameter names and entries correspond to 
%             respective parameter values
% ng - Various scalar model quantities such that
% ng(1) - Number of genes in model trnmodel
% ng(2) - Number of Proteins
% ng(3) - Number of intracellular metabolites
% ng(4) - Number of extracellular metabolites
% solverP - Parameters for SUNDIALS ODE solver
% initSolution - Structure of initial Solutions for the given inputs if 
%                available 
% varname - Names of genes, proteins or metabolites whose solution is to be
%           plotted
% saveData - Struture containing information needed to store text files of
%            data and png files of figures     
%**************************************************************************
%Notes:
%October 2013 
%Added capability to spit out finall steady state values as well instead of
%calling getssval in the function routine
%January 11th 2014
%**************************************************************************
function [initval,Solution,status,finalSS,flag,indx,...
          hallfig,defparval,exptnum] = trnperturbation(model,FBAmodel,batch,...
          defparval,ng,solverP,initSolution,varname,saveData,exptnum)
%The perturbed Output is Passed as an output argument
initflag = 0;
fflag = 1;
if nargin < 9
    exptnum = 0;
end
if nargin < 8
    saveData = struct([]);
end
if nargin < 7
    varname = {};  
    fflag = 0;
end
if nargin < 6
    initSolution = {};
elseif ~isemptyr(initSolution)
    initflag = 1;
end
%ODE Solver Parameters Structure
if nargin < 5
    solverP = struct();
    solverP.tmax = 10;%h
    solverP.tpertmax = 10;%h
    solverP.scalAbsTol = 1e-7;
    solverP.RelTol = 1e-8;
    solverP.MaxIter = 1000;
    solverP.MaxDataPoints = 200;
end
solverP.tmax = batch.tmax;
solverP.tpmax = batch.tpmax;
type = 1; %Type of solution to use for perturbation calculations
%% %SS Solution Initialization
model.nvar = sum(ng)+1;%%[mRNA,Protein,Metabolite]'
if ~initflag %initSolution is empty
    flag = 0;%Initial Status = Not Solved initial SS problem
    fprintf('Evaluating Initial SS Solution \n'); 
    initval = zeros(model.nvar,1);
%     initval(end) = 1.31;%gDCW Non-zero Biomass
    %mRNA umole
    initval(1:ng(1)) = model.SSmRNA;
    %Protein umole
    initval(ng(1)+1:ng(1)+ng(2)) = model.SSreg(1:ng(2));
    %Common Protein-Metabolite
    initval(ng(1)+model.PMind_R) = model.SSreg(ng(2)+1:ng(2)+ng(3));
    initval(ng(1)+model.rnap_ind) = 0;    
    %Intracellular Metabolite mmole
    initval(ng(1)+ng(2)+ng(3)+2:sum(ng)-ng(5)-ng(6)+1) =...
    model.SSreg(ng(2)+ng(3)+2:ng(2)+1+ng(3)+ng(4));
    %Specific Initial Values
    lacind = find(strcmpi('lac[c]',model.Regulators));
    initval(lacind+ng(1)) = 1e-10;
    %Extracellular Metabolite mmole
    model.ext_MC = assign_extconc(batch.init{1},batch.init{2},model);
    initval(sum(ng)-ng(6)+2:sum(ng)+1) = model.ext_MC;
    
    svdata.t = zeros(length(initval),1);
    svdata.y = initval';
    if isfield(saveData,'filename') && ~isempty(saveData.filename)
        savefile(svdata,'initval',saveData);
    else
        fprintf('Not saving file. No file/directory given\n');
    end
    if fflag
        hsfig = figure('Units','pixels',...
                       'Position',[0 50 750 500],...
                       'Visible','off',...
                       'Name','Initial Steady State');  
    else
        hsfig = 0;
    end        
else%Initial Solution is available
    initval = initSolution.y(:,end);
    if isfield(initSolution,'pertbind')
        pardata = struct();
        pardata.pbind = initSolution.pertbind;
    end
    flag = 1;
    if fflag
        hsfig = figure('Units','pixels',...
                       'Position',[0 50 750 500],...
                       'Visible','off',...
                       'Name','Experimental Curve');   
    else
        hsfig = 0;
    end
end
model.Yref = ones(model.nvar,1);
fprintf('Input Condition(s) \n');
fprintf('Input Metabolites\t\t\tConcentrations\n');
for i =1:length(batch.init{1})
    fprintf('%s\t\t\t%3.4g\n',batch.init{1}{i},batch.init{2}(i));
end
fprintf('Evaluation Horizon:%d hours \n',batch.tmax);

pardata.par = model.allpar;
pardata.ng = ng;%
pardata.nvar = model.nvar;
pardata.pdecay = defparval.pdecay;
pardata.mdecay = defparval.mdecay;
pardata.rephill = 3;%defparval.rephill;
pardata.ext_MC = model.ext_MC;
pardata.flux = model.Vss;
pardata.kcat = model.kcat;
pardata.gmax = model.gmax;
pardata.vuptake = model.Vuptake;
pardata.vefflux = model.Vefflux;
pardata.Yref = model.Yref;
%% %Solve ODE System using SUNDIALS
[initSS,status,indx,nfcall,totaltfinish] =...
solveODE(initval,pardata,model,FBAmodel,solverP);
Solution.initSS = initSS;   
if status < 0
    indx = 0;
    finalSS = struct([]);
    hallfig = [];
    return
end
%Obtaining final SS values for initial problem 
f_flag = 0;
for i = 1:size(Solution.initSS.y,2)
    if any(Solution.initSS.y(:,i)<0)
        if i > 1
            finalSS.y = Solution.initSS.y(:,i-1);
            finalSS.t = Solution.initSS.t(i-1);
            f_flag = 1;  
            break
        else
            finalSS.y = Solution.initSS.y(:,i);
            finalSS.t = Solution.initSS.t(i);
        end
        f_flag = 1;    
    end
end
if ~f_flag
    finalSS.t = Solution.initSS.t(end);
    finalSS.y = Solution.initSS.y(:,end); 
end
%% %Initial SS Solution Figure
if hsfig     
    FProperty = struct();        
    FProperty.NumberTitle = 'off';
    FProperty.Color = [0.8 0.8 0.8];           
    [hfig] = dynamicplot(model,ng,varname,initSS,hsfig,[],[],1);        
    set(hsfig,FProperty);
end  
hallfig(1) = 0;
flag = 1; %Status = Solved Initial SS problem 
%Print basic diagnostic
%Advanced disgnostic printed by CVodeMonitor
fprintf('Total Number of Calls to Material Balance Function = %d \n',nfcall);
fprintf('Time for Model Integration using SUNDIALS = %fs \n',totaltfinish);
pertbflag = 0;
%% %Perturbation to Inputs
if pertbflag %If perturbations are defined
    %if length(pertb.metabName) == length(pertb.metabConc)
        npertb = length(batch.pertb{1});
        if flag  %If Initial SS problem has been solved
            for ipertb = 1:npertb                
                fprintf('Performing perturbation #%d out of %d\n',ipertb,npertb);               
                if type == 1
                    %Start Perturbation w/ original initial SS solution
                    %tpertmax = 100000;
                    [pertbSS,pfflag] = perturbationcalc(ipertb,initSS,...
                                       initInConc,solverP,tpertmax);
                    %pfflag = 0;                    
                    %store Solution to Perturbation
                    Solution.(sprintf('pertb%d',ipertb)) = pertbSS;
                    svdata = [pertbSS.t,pertbSS.y'];
                    %pfflag = 0;                    
                    %Obtaining final SS values for perturbation problem                   
                    finalTdata(:,ipertb+1) =...
                        Solution.(sprintf('pertb%d',ipertb)).t(end);
                    finalSS(:,ipertb+1) =...
                        Solution.(sprintf('pertb%d',ipertb)).y(:,end);                 
                    %Perturbed SS Solution Figure                    
                    if pfflag
                        hpfig = figure('Units','pixels',...
                                       'Position',[0 50 2050 750],...
                                       'Visible','off');
                        
                        AProperty.XMinorTick = 'on';
                        AProperty.XTick = 0:1000:tpertmax;
                        AProperty.YMinorTick = 'on';
                        AProperty.TickLength = [0.02,0.02];
                        AProperty.XColor = [.2 .2 .2];
                        AProperty.YColor = [.2 .2 .2];
                        AProperty.XLim = [0 tpertmax];
                        AProperty.FontName = 'Courier';
                        AProperty.FontSize = 12;
                        
                        LgProperty.LineWidth = 1.5;
                        LgProperty.Color = [0 .5 0];                        
                        LpProperty.LineWidth = 1.5;
                        LpProperty.Color = [.8 .4 0];
                        
                        FProperty = struct();
                        figname =...
                        sprintf('#%dPerturbation Response Curves',ipertb);
                        FProperty.Name = figname;
                        FProperty.NumberTitle = 'off';
                        FProperty.Color = [0.2 0.2 0.2];
                        %FProperty.FontSize = 12;
                        %FProperty.FontName = 'Courier';
                        
                        %dbstop in dynamicplot
                        var_name = sprintf('pertb%d',ipertb);
                        hpsubfig = dynamicplot(model,genename,pname,...
                                   Solution.(var_name),hpfig,AProperty,...
                                   LgProperty,LpProperty);                       
                        
                        %hallfig(ipertb+1) = hpfig;
                        
                        set(hpfig,FProperty);
                        
                    else
                        hpfig = 0;
                        
                    end
                    
                    datafname =...
                        sprintf('%s_pertb%d.txt',saveData.filename,ipertb);
                    figfname =...
                        sprintf('%s_pertb%d',saveData.filename,ipertb);                    
                    %savefile(data,datafname,saveData,hfig,figfname)
                    hallfig(ipertb+1) =...
                        savefile(svdata,datafname,saveData,hpfig,figfname);              
                elseif type == 2
                    %Start perturbation with different initial SS solution
                end
                
                %==========================================================
                %Perturbed SS Solution Figure
                %==========================================================           
                
                %Proceed to next perturbation or end
            end         
        else
            %Solve Initial problem first and set flag = 1
            fprintf('UnSolved Initial Problem \n No Initial SS to work with \n');
            return
%             [initSS,status,indx,nfcall,totaltfinish,itime] =...
%                 test3(trnmodel,defparval,ngene,tnreg,ngap,initInConc,tmax);            
        end
else
    fprintf('No Perturbations\n');
%     fprintf('Size Mismatch \n');
    return
    %end 
end
%% Perturbation Function
function [pertbSS,pflag,InConc,nfcall,tfinish] =...
        perturbationcalc(ipertb,initSS,initInConc,solverP,tpertmax)

    if length(batch.pertb{1}{ipertb}) ==...
                length(batch.pertb{2}{ipertb})                     

       %Change to reflect just perturbed metabolites only

       unpertMetab = setdiff(model.Metabolites,batch.pertb{1}{ipertb});
       nmetab = length(model.Metabolites);
       npertmetab = length(batch.pertb{1}{ipertb});
       mflag = 1;

       if length(unpertMetab) ~= nmetab-npertmetab
           fprintf('\n Perturbation Metabolite Mismatch\n');
           fprintf('Cannot find all metaboites in pertb.metabName\n');
           fprintf('Note:Metabolite names are case sensitive\n');                           

           mflag = 0;                   
       end               

       if mflag %If perturbation metabolite names are a match 
%            InConc = inputmetabconc(pertb.metabName{ipertb},...
%                                 pertb.metabConc{ipertb},model);
           imetab = 1;

           while imetab <= nmetab - npertmetab                     

               mbx = strcmp(unpertMetab{imetab},model.Metabolites);
               InConc(mbx) = initInConc(mbx);
               imetab = imetab + 1;
           end               

           fprintf('Total number of perturbed metabolites: %d\n',...
                    length(batch.pertb{1}{ipertb})); 
           fprintf('Perturbed metabolite(s):\n');

           kmetab = 1;
           while kmetab <= length(batch.pertb{1}{ipertb}) 
               fprintf('(%d) %s \n',kmetab,batch.pertb{1}{ipertb}{kmetab});             
               kmetab = kmetab + 1;
           end
           
           fprintf('Integration Horizon:%d\n',tpertmax);
           solverP.tmax = tpertmax;
           %Function Call to integrate Model for each perturbation
           [pertbSS,status,indx,nfcall,tfinish] =...
           solveODE(model,defparval,ngene,tnreg,ngap,InConc,solverP,...
                    initSS.y(:,end));

           flag = 2; %Status = Solved perturbation problem
           pflag = 1;
           %Print basic diagnostic
           %Advanced disgnostic printed by CVodeMonitor
           fprintf('Calls to Material Balance Function = %d\n',nfcall);
           fprintf('Time for Model Integration using SUNDIALS = %f\n\n',tfinish);
       else
           pertbSS.t = [];
           pertbSS.y = [];
           pertbSS.ys = [];
           fprintf('\nTerminating Pertubation Call\n'); 
           pflag = 0;
           flag = 3; %Status = Perturbation problem not solved
                     %Initial solution provided
           return
       end

    else
        fprintf('Mismatch of Concentration & Name Cell Array Sizes of  \n');

        pertbSS.t = [];
        pertbSS.y = [];
        pertbSS.ys = [];
        flag = 1; %Status = Solved initial SS problem
        return
    end              
end
                                    
end



