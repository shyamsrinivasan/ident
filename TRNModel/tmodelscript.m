%Reference for new data file
%2002 Oshima, Molecular Microbiology - Paper on knockouts of TCSs
%2003 Kao & Liao, PNAS, Paper on NCA
%2005 Yang & Liao, Metabolic Engineering, Paper on NCA

clc
clear all
%==========================================================================
%Create model from text file
%==========================================================================
%Two file names:
%file name 1 : Text file for regulatory interactions of all genes in the
%network
%file name 2: Text file regulatory interactions of genes correpsonding to
%transcription factors

flag1 = 0;
if flag1 
    cd('C:\Users\shyam\SkyDrive\Documents\Courses\Modeling Project'); 
    addpath(genpath('C:\Users\shyam\Documents\MATLAB\sundialsTB'));
    addpath(genpath('C:\Users\shyam\SkyDrive\Documents\Courses\Modeling Project\TRN Model'));
    addpath(genpath('C:\Users\shyam\SkyDrive\Documents\Courses\Modeling Project\Other Matlab Central Files'));
    
    rxfname = 'C:\Users\shyam\SkyDrive\Documents\Courses\Modeling Project\TRN Model\TRNtest2.txt';
    regfname = 'C:\Users\shyam\SkyDrive\Documents\Courses\Modeling Project\TRN Model\TRNtest3.txt';
    
    fname1 = 'C:\Users\shyam\SkyDrive\Documents\Courses\Modeling Project\TRN Model\List of Experiments.txt';
    
    saveData.dirname = 'C:\Users\shyam\SkyDrive\Documents\Courses\Modeling Project\Results\TRN Model version 2\New Results';
elseif flag1 == 0
    cd('C:\Users\shyam\Documents\Courses\CHE 1125 Project'); 
    addpath(genpath('C:\Users\shyam\Documents\MATLAB\sundialsTB'));
    addpath(genpath('C:\Users\shyam\Documents\Courses\CHE 1125 Project\TRN Model'));
    addpath(genpath('C:\Users\shyam\Documents\Courses\CHE 1125 Project\Other Matlab Central Files'));
    
    rxfname = 'C:\Users\shyam\Documents\Courses\CHE 1125 Project\TRN Model\TRNtest2.txt';
    regfname = 'C:\Users\shyam\Documents\Courses\CHE 1125 Project\TRN Model\TRNtest3.txt';
    
    fname1 = 'C:\Users\shyam\Documents\Courses\CHE 1125 Project\TRN Model\List of Experiments.txt';
    
    saveData.dirname = 'C:\Users\shyam\Documents\Courses\CHE 1125 Project\Results\TRN Model version 2\New Results';
end

% rxfname = 'C:\Users\shyam\SkyDrive\Documents\Courses\Modeling Project\TRN Model\TRNtest2.txt';
% regfname = 'C:\Users\shyam\SkyDrive\Documents\Courses\Modeling Project\TRN Model\TRNtest3.txt';

regtype = {'(+)';'(-)';'(+/-)'};

%Parameter Values superseded by values given for defparval structure

%=========================
%Define default parameters
%=========================
defparval.accoeff = 1e-9;
defparval.repcoeff = 1e-7;

defparval.rephill = 2;
defparval.dualcoeff = 1e-9;
defparval.brate = 1.66e-6;
defparval.srate = 1.66e-3;%Changed from 1.66e-4;
defparval.trate = 4.94e-6;
defparval.drate = 1.44e-2;%Small oscillations in response seen %Changed from 1.44e-4;
defparval.kmax = 4e-2; %Changed from 5e-3%Changed from 5e-5
defparval.ks = 1e-7;%10 did not provide the expected result. 
%It actually worsened the situation
%Changed from 1e-3;%Changed from 1e-5 Changed from 1e-6
defparval.pdrate = 1.44e-2;%Changed from 1.44e-4 
defparval.lkrate = 1.5e-12; %Transcriptional leakage rate for promoters 
%that are repressed in the absence of an activator

defparval.repcoeff = 1e-7;

% exptnum = 1;

tstart = tic;

%fname1 = 'C:\Users\shyam\Documents\Courses\CHE 1125 Project\TRN Model\List of Experiments.txt';
fileid = fopen(fname1);

if fileid == -1
    fprintf('File %s cannot be opened.', fname1);
    trnmodel = struct([]);
    return;
end

E = textscan(fileid, '%s%s%s%s%s%s%s', 'Delimiter', '\t',...
            'TreatAsEmpty', {'None'}, 'HeaderLines', 1);
fclose(fileid);

%newterms = regexp(xterms{1},'(\w+.?)(\(\W+.?\))+','tokens');
if ~isempty(E)
    nexpt = length(E{1});
else
    nexpt = 0;
end
allSolution = struct();
SSsolution = struct();
for exptnum = 1:1
    fprintf('\nBeginning Simulation #%d of %d\n',exptnum,nexpt);   
    %==========================================================================
    %Initialize Model using specific/default values
    %==========================================================================
    
    [defparval] = ExtractParameters(E{7}(exptnum),defparval);
    %dbstop in Tmodel.m
    [trnmodel,defparval,ngene,ngap,tnreg,regprotein,rgene,rmetab,...
        rprot,defvalueGene,C] = Tmodel(rxfname,regfname,defparval);
    
    GeneStruct = struct();
    GeneStruct.regulators = {{'ArcA';'FNR'};{'FNR';'OxygenRecp'}};%Not of consequence currently
    GeneStruct.coeff = {[1e-7;1e-7];[1e-5;1e-5]};%Not of consequence currently
    GeneStruct.srate = {[1e-3;1e-4];[1e-4;1e-4]};%Not of consequence currently
    %==========================================================================
    %Change individual/all parameters in a given model using changeparameter.m
    %==========================================================================
    [trnmodel] = changeparameter(trnmodel,defparval,regprotein,GeneStruct);
    %==================================================================
    %Blacklisting unwanted transcription factor activities
    %==================================================================    
    %if No regulation has been specified for a regulatory protein,
    %Blacklist correposnding protein from affecting transcription of
    %other genes
    %Added December 05th 2013    
    % blacklist = {'CpxR';'FIS';'IHF';'NarL';'RpoN';'RpoS';'SoxS'};
    blacklist = {};
    if ngap == length(regprotein) && ngap == length(rgene) &&...
            ngap == length(rmetab) && ngap == length(rprot)
        iprot = 1;
        while iprot <= ngap
            if ~isempty(regprotein{iprot}) && ~isempty(rgene{iprot})
                if isempty(rmetab{iprot}) && isempty(rprot{iprot})
                    blacklist = [blacklist;regprotein{iprot}];
                end
            end
            iprot = iprot + 1;
        end
    end
    
    [trnmodel,statusflag] = blacklist_tfs(trnmodel,blacklist,regprotein);
    %==========================================================================
    %Define Initial Conditions/Parameters and Solve ODEs
    %==========================================================================
    %inMetab & inMetabConc are returned as cells
    dbstop in ExtractData.m
    [inMetab,inMetabConc] = ExtractData(E{2}(exptnum),E{3}(exptnum));
    [pertMetab,pertMetabConc] = ExtractData(E{4}(exptnum),E{5}(exptnum));
    
    metabName = inMetab{1};
    metabConc = inMetabConc{1};
    
    pertb.metabName = pertMetab;
    pertb.metabConc = pertMetabConc;     
   
    compos = strfind(E{6}{exptnum},',');
    ngenenames = length(compos) + 1;
    genename = cell(ngenenames,1);
    
    for iname = 1:ngenenames
        if ngenenames > 1 && iname == ngenenames
            genename{iname,1} = E{6}{exptnum}(compos(iname-1)+1:end);
        elseif ngenenames > 1 && iname > 1
            genename{iname,1} = E{6}{exptnum}(compos(iname-1)+1:compos(iname)-1);
        elseif ngenenames > 1 && iname == 1
            genename{iname,1} = E{6}{exptnum}(1:compos(iname)-1);
        else%ngenenames = 1
            genename{iname,1} = E{6}{exptnum}(1:end);
        end
    end 
    
    pname = {};
    %pname = {'OxygenRecp';'GLCxtRecp'};
    
    initSolution = {};
    
    %CVode Solver Parameters
    solverP.tmax = 2500;
    solverP.scalAbsTol = 1e-5;
    solverP.RelTol = 1e-5;
    solverP.MaxIter = 1000;    
    solverP.MaxDataPoints = 200;    
    tpertmax = 2500;      
    
    saveData.filename = sprintf('ExptCondition_%d',exptnum);    
    %save.dirname = 'C:\Users\shyam\Documents\Courses\CHE 1125H\Project\Results';
    %saveData.Figfilename = sprintf('initSS_Experiment_%d_%s',exptnum,datestr(date));
      
    trnmodel.allpar = parameter_vector(trnmodel.Coefficient,trnmodel.srate,ngene);
    fields = {'Coefficient','srate'};
    trnmodel = rmfield(trnmodel,fields);
    trnmodel = rmfield(trnmodel,'Regulator');
    
    
    [InConc,Solution,finalSSdata,finalTdata,flag,odestatus,indx,figh] =...
        trnperturbation(trnmodel,defparval,ngene,tnreg,ngap,...
        inMetab{1},inMetabConc{1},solverP,initSolution,pertb,tpertmax,...
        genename,pname,saveData,exptnum);    
    
    
    if exist('Solution','var')
        genename = {'arcA','fnr','sdhA','aceA'};
        pname = {'OxygenRecp'};
        %dbstop in getssval.m
        %dbstop in storeSolution.m
        %[finalSSdata,finalTdata,exptname] = getssval(Solution,ngene,tnreg);
        [allSolution,SSsolution] =...
         storeSolution(Solution,SSsolution,allSolution,trnmodel,ngene,genename,pname);
        %genename = {'fnr'};
        %[allSolution] = storeSolution(Solution,allSolution,trnmodel,ngene,tnreg,genename,pname);
        %storeSolution(Solution,AllDynSolution,trnmodel,gname,pname,exptnum)
    end
    
    %Ploting inidividual SS values for different conditions
%     genename = {'arcA';'fnr'};
%     
%     pname = {};
    %dbstop in plotSSexpression.m
    %plotSSexpression(trnmodel,finalSSdata,genename,pname,ngene,tnreg,saveData);
    
    %PLot dynamic profiles for different perturbations for a given
    %experiment
    
     close all        
    
    fprintf('Completed Simulation #%d of %d\n',exptnum,nexpt);
    
end

tend = toc(tstart);
fprintf('Total Time for Model Simulation %fs\n\n',tend);

%%


%Log transform & Noise addition
Data = Solution.initSS{2};
Options.MaxDataPoints = solverP.MaxDataPoints;
Options.NoiseSD = 5e-1;
[logdata,LogNoisyData] = NoiseTransform(Data,Options);

%Checking for noisy data spread
logfig = plot(Solution.initSS{1}(1:100),LogNoisyData(1,1:100));
% figh = plot(Solution.initSS{1}(1:100),NoisyData(1,1:100))
set(logfig,'LineStyle','none','Marker','.','Color',[.3 .3 .3])
hold on
figh = plot(Solution.initSS{1}(1:100),logdata(1,1:100));
hold off
%%
%======================================================
%Simulating smaller sub-system models of TRN
%======================================================


%Create large model from text file
flag1 = 0;
if flag1 
    cd('C:\Users\shyam\SkyDrive\Documents\Courses\Modeling Project'); 
    addpath(genpath('C:\Users\shyam\Documents\MATLAB\sundialsTB'));
    addpath(genpath('C:\Users\shyam\SkyDrive\Documents\Courses\Modeling Project\TRN Model'));
    addpath(genpath('C:\Users\shyam\SkyDrive\Documents\Courses\Modeling Project\Other Matlab Central Files'));
    
    rxfname = 'C:\Users\shyam\SkyDrive\Documents\Courses\Modeling Project\TRN Model\TRNtest2.txt';
    regfname = 'C:\Users\shyam\SkyDrive\Documents\Courses\Modeling Project\TRN Model\TRNtest3.txt';
    
    fname1 = 'C:\Users\shyam\SkyDrive\Documents\Courses\Modeling Project\TRN Model\newList of Experiments.txt';
    
    saveData.dirname = 'C:\Users\shyam\SkyDrive\Documents\Courses\Modeling Project\Results\TRN Model version 2\New Results';
elseif flag1 == 0
    cd('C:\Users\shyam\Documents\Courses\CHE 1125 Project'); 
    addpath(genpath('C:\Users\shyam\Documents\MATLAB\sundialsTB'));
    addpath(genpath('C:\Users\shyam\Documents\Courses\CHE 1125 Project\TRN Model'));
    addpath(genpath('C:\Users\shyam\Documents\Courses\CHE 1125 Project\Other Matlab Central Files'));
    
    rxfname = 'C:\Users\shyam\Documents\Courses\CHE 1125 Project\TRN Model\TRNtest2.txt';
    regfname = 'C:\Users\shyam\Documents\Courses\CHE 1125 Project\TRN Model\TRNtest3.txt';
    
    fname1 = 'C:\Users\shyam\Documents\Courses\CHE 1125 Project\TRN Model\newList of Experiments.txt';
    
    saveData.dirname = 'C:\Users\shyam\Documents\Courses\CHE 1125 Project\Results\TRN Model version 2\New Results';
end
%Parameter Values superseded by values given for defparval structure
%=========================
%Define default parameters
%=========================
defparval.accoeff = 1e-9;
defparval.repcoeff = 1e-7;
defparval.rephill = 2;
defparval.dualcoeff = 1e-9;
defparval.brate = 1.66e-6;
defparval.srate = 1.66e-3;%Changed from 1.66e-4;
defparval.trate = 4.94e-6;
defparval.drate = 1.44e-2;%Small oscillations in response seen %Changed from 1.44e-4;
defparval.kmax = 4e-2; %Changed from 5e-3%Changed from 5e-5
defparval.ks = 1e-7;%10 did not provide the expected result. 
%It actually worsened the situation
%Changed from 1e-3;%Changed from 1e-5 Changed from 1e-6
defparval.pdrate = 1.44e-2;%Changed from 1.44e-4 
defparval.lkrate = 1.5e-12; %Transcriptional leakage rate for promoters 
%that are repressed in the absence of an activator
defparval.repcoeff = 1e-7;

% exptnum = 1;

tstart = tic;
fileid = fopen(fname1);

if fileid == -1
    fprintf('File %s cannot be opened.', fname1);
    trnmodel = struct([]);
    return;
end

E = textscan(fileid, '%s%s%s%s%s%s%s', 'Delimiter', '\t',...
            'TreatAsEmpty', {'None'}, 'HeaderLines', 1);
fclose(fileid);

if ~isempty(E)
    nexpt = length(E{1});
else
    nexpt = 0;
end
allSolution = struct();
SSsolution = struct();
for exptnum = 2:2
    fprintf('\nBeginning Simulation #%d of %d\n',exptnum,nexpt);      
    %Initialize Model using specific/default values    
    [defparval] = ExtractParameters(E{7}(exptnum),defparval);    
    [trnmodel,defparval,ngene,ngap,tnreg,regprotein,rgene,rmetab,...
        rprot,defvalueGene,C] = Tmodel(rxfname,regfname,defparval);     
        
    %Change individual/all parameters in a given model 
    GeneStruct = struct();
    [trnmodel] = changeparameter(trnmodel,defparval,regprotein,GeneStruct);
    
    %Blacklisting unwanted transcription factor activities
    blacklist = {};
    if ngap == length(regprotein) && ngap == length(rgene) &&...
            ngap == length(rmetab) && ngap == length(rprot)
        iprot = 1;
        while iprot <= ngap
            if ~isempty(regprotein{iprot}) && ~isempty(rgene{iprot})
                if isempty(rmetab{iprot}) && isempty(rprot{iprot})
                    blacklist = [blacklist;regprotein{iprot}];
                end
            end
            iprot = iprot + 1;
        end
    end    
    [trnmodel,statusflag] = blacklist_tfs(trnmodel,blacklist,regprotein);
    
    %Define Initial Conditions/Parameters and Solve ODEs
    [inMetab,inMetabConc] = ExtractData(E{2}(exptnum),E{3}(exptnum));
    [pertMetab,pertMetabConc] = ExtractData(E{4}(exptnum),E{5}(exptnum));
    
    metabName = inMetab{1};
    metabConc = inMetabConc{1};
    
    pertb.metabName = pertMetab;
    pertb.metabConc = pertMetabConc;  
    
    compos = strfind(E{6}{exptnum},',');
    ngenenames = length(compos) + 1;
    genename = cell(ngenenames,1);
    
    for iname = 1:ngenenames
        if ngenenames > 1 && iname == ngenenames
            genename{iname,1} = E{6}{exptnum}(compos(iname-1)+1:end);
        elseif ngenenames > 1 && iname > 1
            genename{iname,1} = E{6}{exptnum}(compos(iname-1)+1:compos(iname)-1);
        elseif ngenenames > 1 && iname == 1
            genename{iname,1} = E{6}{exptnum}(1:compos(iname)-1);
        else%ngenenames = 1
            genename{iname,1} = E{6}{exptnum}(1:end);
        end
    end     
    pname = {};
    initSolution = {};
    
    %CVode Solver Parameters
    solverP.tmax = 2500;
    solverP.scalAbsTol = 1e-5;
    solverP.RelTol = 1e-5;
    solverP.MaxIter = 1000;    
    solverP.MaxDataPoints = 200;    
    tpertmax = 2500; 
    
    %other parameters
    saveData.filename = sprintf('ExptCondition_%d',exptnum);
    %save.dirname = 'C:\Users\shyam\Documents\Courses\CHE 1125H\Project\Results';
    %saveData.Figfilename = sprintf('initSS_Experiment_%d_%s',exptnum,datestr(date));
    
    trnmodel.allpar = parameter_vector(trnmodel.Coefficient,trnmodel.srate,ngene);
    fields = {'Coefficient','srate'};
    trnmodel = rmfield(trnmodel,fields);
    trnmodel = rmfield(trnmodel,'Regulator');
    
    %generating smaller subsystem model
    gene = {'aceA'};
    [newtrnmodel,regprotein] = separatesubsystems(gene,trnmodel,regprotein);
    ngene = length(newtrnmodel.Gene);
    tnreg = length(newtrnmodel.Protein);
    ngap = length(newtrnmodel.Protein(...
                  cellfun('isempty',regexp(newtrnmodel.Protein,'\w+.?Recp'))));

    [coefficient,srate] = parameter_return(newtrnmodel.allpar,newtrnmodel,ngene,tnreg);
    newtrnmodel.Coefficient = coefficient;
    newtrnmodel.srate = srate;
    blacklist = {'FruR';'IclR';'CRP';'FadR';'exFDPRecp';'exF6PRecp'};
    [newtrnmodel,statusflag] = blacklist_tfs(newtrnmodel,blacklist,regprotein);

    newtrnmodel.allpar = parameter_vector(newtrnmodel.Coefficient,newtrnmodel.srate,ngene);
    fields = {'Coefficient','srate'};
    newtrnmodel = rmfield(newtrnmodel,fields);    
    
    [InConc,Solution,finalSSdata,finalTdata,flag,odestatus,indx,figh] =...
        trnperturbation(newtrnmodel,defparval,ngene,tnreg,ngap,...
        metabName,metabConc,solverP,initSolution,pertb,tpertmax,...
        genename,pname,saveData,exptnum);    
    
    
    if exist('Solution','var')
        genename = {'arcA','fnr','sdhA','aceA'};
        pname = {'OxygenRecp'};
        
        [allSolution,SSsolution] =...
         storeSolution(Solution,SSsolution,allSolution,newtrnmodel,ngene,genename,pname);        
    end   

     close all        
    
    fprintf('Completed Simulation #%d of %d\n',exptnum,nexpt);
    
end

tend = toc(tstart);
fprintf('Total Time for Model Simulation %fs\n\n',tend);

%Generate Pseudo Experimental Data by adding noise to model simulation
%results
Data = Solution.initSS{2};
Options.MaxDataPoints = solverP.MaxDataPoints;
Options.NoiseSD = 5e-1;
[logdata,LogNoisyData,Noise] = NoiseTransform(Data,Options);

% expt.t = Solution.initSS{1};
% expt.y = LogNoisyData;
% expt.error = Noise;

calc.t = Solution.initSS{1};
calc.y = logdata;

%[Output] = LogNormPDF( Values, Means, Variance);

%%
%==========================================================================
%Plotting Dynamic data obtained from storeSolution.m
%==========================================================================
genename = {'arcA','fnr','sdhA'};
pname = {'OxygenRecp'};

%dbstop in plotSSexpression.m
saveData.filename = 'Comparisons';
plotSSexpression(trnmodel,[],genename,pname,ngene,tnreg,saveData,allSolution,[15 20]);


%%
%Extract Sensitivity Information from Solution.x
%YS = Solution.initSS{3};

%Include This in trnperturbation.m

Ns = 2;
[N,ncolumns] = size(Solution.initSS{3});
ntimepoints = ncolumns/Ns;
%YSnew = zeros(ntimepoints,N,Ns);
Sens = cell(Ns,1);
%ipar*icolumn - (npar-ipar)
for ipar = 1:Ns
    Sens{ipar} = zeros(N,ntimepoints);
    for itime = 1:ntimepoints
        Sens{ipar}(:,itime) = Solution.initSS{3}(:,Ns*itime-(Ns-ipar));
    end
end

YS = zeros(N,ntimepoints,Ns);
%YS(:,:,1) = Sens{1};
for idim = 1:Ns
    YS(:,:,idim) = Sens{idim};
end
%%
%==========================================================================
%Function for Likelihood and Log-likelihood estimation
%Copied from LogNormPDF of the Bayesian ED Technique
%Written as LogNormPDF.m
%April 04 2014
%==========================================================================
%function for determination of log-likelihood
%Liklehood
%L(Y|Theta) = PI(1/sqrt(2pi)sigma)*EXP(-(Data(i,j) - Model(i,j))/2sigma)
%Variables for likelihood for gaussian variable
%Expt value - y(t)
%Model value or mean - x(t)
%Variance - sigma
%NoisyData = NoisyData + randn(NumOfSpecies, D).*NoiseSD;

%Output = sum( -ones(D,1)*(0.5*log(2*pi*Variance)) - ((Values-Means).^2)./(2*(ones(D,1)*Variance))
%Values - Estimates
%Means - Data

variance = ones(1,1)*5e-1;
Data = LogNoisyData(1,:);
ModelOutput = Solution.initSS{2}(1,:);
% ES2 = ((ModelOutput - Data).^2)';
% DN = 2*ones(200,1)*variance;
% Coeff = ones(200,1)*0.5*log(2*pi*variance);
% LH = sum(-Coeff - ES2./DN);
%LH = sum(-ones(200,1)*(0.5*log(2*pi*variance))-((ModelOutput - Data).^2)./(2*(ones(200,1)*variance)));

LH = LogNormPDF(ModelOutput,Data,variance);


%%

%==========================================================================
%Solving for problems that already have a initial Solution vector available
%==========================================================================

% if exist('Solution','var')
%     if isfield(Solution,'pertb1_SS')
%         if ~isemptyr(Solution.pertb1_SS)
%             %Initial Solution is present in workspace and can be
%             %used
%             %dbstop in trnperturbation
%             %Using perturbation solution from anaerobic to aerobic
%             [InConc,Solution,flag,odestatus,indx,figh,initdefparval] =...
%             trnperturbation(newtrnmodel,defparval,ngene,tnreg,ngap,...
%             metabName,metabConc,tmax,Solution.pertb1_SS,pertb,tpertmax,...
%             genename,pname);
%         end
%     end
% else
%     %No Initial SS available solve for initial SS
%     fprintf('No Initial Solution. Abort!\n');
% %     initSolution = {};
% %     [InConc,Solution,flag,odestatus,indx,figh] =...
% %     trnperturbation(newtrnmodel,defparval,ngene,tnreg,ngap,...
% %     metabName,metabConc,tmax,initSolution,pertb,tpertmax,...
% %     genename,pname);       
% end


% defparval.srate = 1e-4/(10^jrate);
%     end
% defparval.repcoeff = 1e-9*(10^irate);
% end

%%
dbfname = 'C:\Users\shyam\Documents\Courses\CHE 1125H\Project\TRNtest2.xlsx';
sheetnum = 'Sheet4';
wflag = 1;
Regulators = C{2};
pflag = 0;
%dbstop in srateinit.m
[strates,btrate,pflag] = srateinit(trnmodel,Regulators,defparval,pflag,wflag,dbfname,sheetnum);


%%
%Setting axis labels using object handles 

set(get(h(2),'XLabel'),'FontName','Courier')
set(gca,'Title',text('String','New Title','Color','r'))

%%
%Consistency Check code
clc
coeffsize = zeros(ngene,1);
for igene = 1: ngene
    if length(find(trnmodel.RS(igene,:)))==length(bindaff{igene})
        coeffsize(igene) = 1;
    end    
end   
    
    
% check = zeros(ngene,1);
% for igene = 1:ngene
%     if isequal(find(trnmodel.Order(igene,:)),find(trnmodel.Coefficient(igene,:)))
%         %check(igene) = 1;
%     end
%     if length(find(trnmodel.RS(igene,:))) == 1
%         if trnmodel.Order(igene,:) ~= 1
%             %check(igene) = 1;
%         end
%     end
%     if isequal(find(trnmodel.RS(igene,:)),find(trnmodel.srate(igene,:)))
%         check(igene) = 1;
%     end
% end


%%
%Miscellaneous RegExp function calls

%[proterms,xterms] = regexp(terms{iterm},'(\w+.?)(\W+.?)\[(\w*.?)\]+','tokens','split');
%[nproterms] = regexp(xterms{end},'(\w+.?)(\W+.?)','tokens');

clc
var = {'ArcA-|FruR-|IclR-|IHF+'};
[terms,xterms] = regexp(var{1},'(\w+.?)(\W+?)(\W+?)+','tokens','split');
nterms = regexp(xterms{end},'(\w+.?)(\W+?)+','tokens');
terms = [terms,nterms];

clc
%[Srate] = splitReg(trnmodel,ngene);



%%

   