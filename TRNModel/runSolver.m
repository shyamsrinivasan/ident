% Wrapper to call solveODE or solveDAE
function [initval,Solution,status,finalSS,flag,indx,...
          hallfig,defparval,exptnum] = runSolver(model,FBAmodel,batch,...
          defparval,ng,solverP,initSolution,varname,saveData,exptnum)
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
%% %SS Solution Initialization
model.nvar = sum(ng)+1;%%[mRNA,Protein,Metabolite]'
if ~initflag %initSolution is empty
    flag = 0;%Initial Status = Not Solved initial SS problem
    fprintf('Evaluating Initial SS Solution \n'); 
    initval = zeros(model.nvar,1);
    initres = zeros(model.nvar,1);
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

%Variable type
%0 - Algebraic, 1 - Differential
pardata.varind = ones(model.nvar,1);
algind = ~cellfun('isempty',strfind(model.Regulators,'-')) |...
         ~cellfun('isempty',strfind(model.Regulators,'*')) |...
         ~cellfun('isempty',regexp(model.Regulators,'^(?:P_)\w'));
pardata.varind([false(ng(1),1);algind]) = 0;     
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
%% %Solve DAE System using SUNDIALS
if ~all(pardata.varind)
    [initSS,status,indx,nfcall,totaltfinish] =...
    solveDAE(initval,initres,pardata,model,FBAmodel,solverP);  
else
    %Solve ODE System using SUNDIALS
    [initSS,status,indx,nfcall,totaltfinish] =...
    solveODE(initval,pardata,model,FBAmodel,solverP);   
end
Solution.initSS = initSS;     
if status < 0
    indx = 0;
    finalSS = struct([]);
    hallfig = [];
    return
end
%Obtaining final SS values for initial problem 
f_flag = 0;
i = size(Solution.initSS.y,2);
while ~f_flag
% for i = size(Solution.initSS.y,2):-1:1
    if any(Solution.initSS.y(1:ng(1),i)<-solverP.RabsTol) ||...
       any(Solution.initSS.y(ng(1)+1:ng(1)+ng(2)+ng(3),i)<-solverP.PabsTol)||...
       any(Solution.initSS.y(ng(1)+ng(2)+ng(3)+2:ng(1)+ng(2)+ng(3)+ng(4),i)<-solverP.MabsTol)
        finalSS.y = Solution.initSS.y(:,i-1);
        finalSS.t = Solution.initSS.t(i-1);
        f_flag = 1;  
    else
        finalSS.y = Solution.initSS.y(:,i);
        finalSS.t = Solution.initSS.t(i);
        f_flag = 1;
    end
    i = i-1;
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
return
