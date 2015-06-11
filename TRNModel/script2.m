%Compare different cases using plots
fname1 = 'C:\Users\shyam\Documents\Courses\CHE 1125 Project\TRN Model\ListExp.txt';
tstart = tic;
fileid = fopen(fname1);
if fileid == -1
    fprintf('File %s cannot be opened.', fname1);
    trnmodel = struct([]);
    return;
end

E = textscan(fileid, '%d%s%s%s%f', 'Delimiter', '\t',...
            'TreatAsEmpty', {'None'}, 'HeaderLines', 1);
fclose(fileid);
floc = 'C:\Users\shyam\Documents\Courses\CHE 1125 Project\TRN Model';
if ~isempty(E)    
    allSolution = struct();
    nexpt = length(E{1});
    for iexpt = 1:nexpt
        %Load initial model
        fprintf('Experiment %d',E{1}(iexpt));
        fname = [floc sprintf('\\%s.mat',E{2}{iexpt})];
        load(fname);
        %or build model
%         [trnmodel,FBAmodel,defparval] = Tmodel(trnfname,regfname,FBAmodel,variable); 
        %Model Initialization
        batch = struct();
        batch.init{1} = {'M1xt';'M2xt'};
        batch.init{2} = [1e3;10];%mmoles
        batch.tmax = 300000;%s
        batch.tpmax = 10;%h  
        %Assign plot variables
        E{3}{iexpt} = strtrim(strrep(E{3}{iexpt},'"',''));
        [varname] = ExtractData(2,E{3}{iexpt});
        %assign Parameters
        [trnmodel] = ExtractData(3,trnmodel,E{4}{iexpt},E{5}(iexpt));
        %ODE solver parameters
        initSolution = {};
        solverP.RabsTol = 1e-6;
        solverP.PabsTol = 1e-6;
        solverP.MabsTol = 1e-6;
        solverP.RelTol = 1e-3;
        solverP.MaxIter = 1000;    
        solverP.MaxDataPoints = 200;         
        saveData.filename = '';%sprintf('ExptCondition_%d',exptnum);
        %Solve Model ODE 
        [initval,Solution,finalSSdata,finalTdata,flag,odestatus,indx,figh] =...
        trnperturbation(trnmodel,FBAmodel,batch,defparval,ng,...
                        solverP,initSolution,...
                        varname,saveData);    
        allSolution.(sprintf('sim_%d',iexpt)) = Solution.initSS;       
    
         close all        

        fprintf('Completed Simulation #%d of %d\n',iexpt,nexpt);
    
    end
end
%Plot Comparisons
plotComparison(trnmodel,nexpt,ng,varname,allSolution);

tend = toc(tstart);
fprintf('Total Time for Model Simulation %fs\n\n',tend);