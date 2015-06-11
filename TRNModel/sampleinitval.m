% function [allSolution,allfinalSS,trnmodel,plotData,h_subfig] =...
%          plotComparison(fname,ng,varname,varargin)
% Calculate and plot comparison for different simulation environments 
% steady state and time course profiles of the integrated network
function [allSolution,allfinalSS,trnmodel,conc,flux,petconc] =...
         sampleinitval(fname,ng,varname,varargin)
     
allSolution = varargin{1};
trnmodel = varargin{2};
parName = varargin{3};
parScale = varargin{4};
flag = 0;
if isempty(allSolution)
    flag = 1;
end
fileid = fopen(fname);
if fileid == -1
    fprintf('File %s cannot be opened.', fname);
    trnmodel = struct([]);
    return;
end

E = textscan(fileid, '%d%s%s%s%s%s%s', 'Delimiter', '\t',...
            'TreatAsEmpty', {'None'}, 'HeaderLines', 1);
fclose(fileid);
floc = 'C:\Users\shyam\Documents\Courses\CHE 1125 Project\TRN Model';
%Colors
load('C:\Users\shyam\Documents\Courses\CHE 1125 Project\Colors.mat');

if flag
    if ~isempty(E)    
        allSolution = struct();
        nexpt = length(E{1});
        par_ = cell(nexpt,2);  
        initSolution = struct([]);
        for iexpt = 1:1
            %Load initial model
            fprintf('Experiment %d\n',E{1}(iexpt));
            fname = [floc sprintf('\\%s.mat',E{2}{iexpt})];
            load(fname);
            fprintf('Model Loaded: %s\n',E{2}{iexpt});  
            %Change Parameters for Sensitivity Analysis
            nsampl = 1;
%             parVal = MCsample_parameters(parName,parScale,nsampl);           
            %for each parameter sample simulate model + MC over initial
            %conditions            
            for jsampl = 1:nsampl
                 %Assign parameters
%                  [trnmodel] = sensAssign(trnmodel,parName,parScale,parVal,jsampl);                 
                %or build model
        %         [trnmodel,FBAmodel,defparval] = Tmodel(trnfname,regfname,FBAmodel,variable); 
                %Model Initialization
                sim_name = sprintf('sim_%d',iexpt);
                %Assign plot variables
                E{3}{iexpt} = strtrim(strrep(E{3}{iexpt},'"',''));
                E{4}{iexpt} = strtrim(strrep(E{4}{iexpt},'"',''));
                E{5}{iexpt} = strtrim(strrep(E{5}{iexpt},'"',''));
                [varname] = ExtractData(2,E{3}{iexpt});                
                if ~isempty(E{4}{iexpt}) && isempty(E{5}{iexpt})
                    %sample initial values experiments
                    ninitsample = ExtractData(4,E{4}{iexpt});
                end
                
                if isempty(initSolution)
                    [trnmodel,batch,solverP,saveData] = initializeModel(trnmodel);
                    fprintf('Cell Growth Rate = %1.2g h-1\n',trnmodel.gmax);
                    %calculate solution for model at gmax and set initSolution 
%                     [~,Solution,~,finalSS] =...
%                     runSolver(trnmodel,FBAmodel,batch,defparval,ng,...
%                               solverP,initSolution,varname,saveData); 
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
                
                ninitsample = 10000;
                conc = zeros(sum(ng)+1,ninitsample+1);
                petconc = zeros(sum(ng)+1,ninitsample);
                flux = zeros(trnmodel.nt_rxn,ninitsample+1);
                y0 = allfinalSS.init.y;
                
                %Perturb initial values
                %Y0 range                
                lb = [1e-4 0.1];
                ub = [1000 100];
                pvar = {'P6','P'};
                [y0new,pbind,savesimd] =...
                MCsample_initval(trnmodel,ng,pvar,y0,lb,ub,ninitsample);                
                
                for isamp = 1:ninitsample
                    fprintf('\nStarting Sample %d of %d \n',isamp,ninitsample);                    
                    ksampl = sprintf('sample_%d',isamp);
                    [trnmodel,batch,solverP,saveData] = initializeModel(trnmodel);
                    fprintf('Cell Growth Rate = %1.2g h-1\n',trnmodel.gmax);
                    saveData.filename = ksampl;
                    %Choose a sampled value
                    petconc(:,isamp) = y0new(:,isamp);
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
                    pertbinitval.y = y0new(:,isamp);
                    pertbinitval.pertbind = pbind;                       
                    trnmodel.ext_MC = [];
                    [~,Solution,~,finalSS] =...
                    runSolver(trnmodel,FBAmodel,batch,defparval,ng,...
                              solverP,pertbinitval,varname,saveData); 
%                     [~,Solution,~,finalSS] =...
%                     trnperturbation(trnmodel,FBAmodel,batch,defparval,ng,...
%                                     solverP,pertbinitval,...
%                                     varname,saveData); 
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
                    fprintf('Completed Sample #%d of %d\n',isamp,ninitsample);
                    close all  
                end  
                fprintf('Completed Parameter Sample #%d of %d\n',jsampl,nsampl);
            end
            fprintf('Completed Simulation #%d of %d\n',iexpt,nexpt);
        end
    end
end

conc(:,1) = allfinalSS.init.y;
flux(:,1) = intflux;
%All figures and Data Collection
printvar = {'P1','P2','P3','P4','P5','P6','P7','P8'};
%write data to excel file
% selectData(allSolution,trnmodel,ng);
%Plot Fluxes & Bin and plot Flux Distribution
plotflux(trnmodel,flux,printvar);
% plotfluxdistr(trnmodel,flux,printvar);
%Plot Steady State Concentrations
plotSSexpression(trnmodel,ng,conc,petconc,varname)
%Plot Solution Curves
% [plotData,h_subfig] =...
% printResults(trnmodel,ng,allSolution,conc,petconc,flux,printvar,varname);
return

function [plotData,hb_subfig] =...
         bifurcationPlot(ColorSpec,par_,nexpt,trnmodel,ng,varname,allfinalSS)
plotData.val = zeros(nexpt,0);
for iexpt = 1:nexpt
    %Select model
    sim_name = sprintf('sim_%d',iexpt);
    %Select Data 
    data.x = str2double(par_{iexpt,1}{2});
    data.y = allfinalSS.(sim_name).y;    
    data.xlabel = par_{iexpt,1}{1};
    %Plot solution curve
    LineP.Marker = 'o';
    LineP.MarkerEdgeColor = 'none';
    LineP.MarkerFaceColor = ColorSpec{iexpt};
    LineP.Displayname = sprintf('%s=%s',par_{iexpt,1}{1},par_{iexpt,1}{2});
    if isempty(findobj('type','figure','Name','Bifurcation Diagram'))
        hb_fig = figure('Name','Bifurcation Diagram','Color',[1 1 1]);
    else
        hb_fig = findobj('type','figure','Name','Bifurcation Diagram');
        figure(hb_fig);
    end
    if exist('hb_subfig','var')        
        [hb_subfig,new_data] =...
        steadystateplot(trnmodel,ng,varname,data,hb_fig,LineP,hb_subfig);
    else
        [hb_subfig,new_data] =...
        steadystateplot(trnmodel,ng,varname,data,hb_fig,LineP);
    end  
    plotData.val(iexpt,1) = new_data.x;
    plotData.val(iexpt,2:length(new_data.y)+1) = new_data.y';    
end
plotData.label = [data.xlabel,new_data.labels'];
setProperties(hb_fig,hb_subfig,[],plotData);
return

function ColorSpec = setLineP(cl_name)
n_color = length(cl_name);
ColorSpec = cell(n_color,1);
for icolor = 1:n_color
%     fac = icolor/(1+icolor);
%     ColorSpec{icolor} = [1/2*fac 2/3*fac 1/3*fac];  
    ColorSpec{icolor} = rgb(cl_name{icolor});
end
    