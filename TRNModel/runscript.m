%==========================
%PL-based ED
%==========================
%Experimental design using PL
%Select aceA-ArcA subsystem model.RS(1,1)
%clear newtrnmodel

%Use pre-defined parameter values (from expt list file)
%for generating pseudo experimental data
%generate a different parameter set for validating ED

%start all the way from reading data and extracting it from the file

% defparval.repcoeff = 1e-5;
% GeneStruct = struct();
% [newtrnmodel] = changeparameter(newtrnmodel,defparval,regprotein,GeneStruct);

%All other parameters are pre-specified from the earlier simulation






% yexpt = expt.y;
% ycalc = calc.y;
% yerror = expt.error;

%Use initial calculated dataset to calculate initial objective function
chi2gene = zeros(ngene,1);
for igene = 1:ngene      
    chi2scal = zeros(length(expt.t),1);
    yexpt = expt.y(igene,:);
    ycalc = calc.y(igene,:);
    yerror = expt.error(igene,:);
    
    iexptime = 1; 
    icalctime = 1;
    while iexptime <= length(expt.t) && icalctime <= length(calc.t)                
        chi2scal(iexptime) = ((yexpt(iexptime)-ycalc(icalctime))/yerror(iexptime))^2;                
        icalctime = icalctime + 1;                
        iexptime = iexptime + 1;
    end
    %Use chi2gene if a local sensitivity measure is required to estimate
    %parameters
    chi2gene(igene) = sum(chi2scal);
end
objval = sum(chi2gene);

%InitConc should also be passed as an input argument to
[psple,objple] = estimatepl(expt,p,obj);



% 
% p = newtrnmodel.allpar;
% 
% %Initial Objective calculation
% %model evaluation
% newtrnmodel.allpar = p;
% 
% [calc,status] =...
% solveODE(newtrnmodel,defparval,ngene,tnreg,ngap,initInConc,solverP);
% 
% %Calculate objective function (a least square estimate)
% 
% chi2gene = zeros(ngene,1);
% 
% exptime = expt.t;
% yexp = expt.y(1:ngene,:);
% yerror = expt.error(1:ngene,:);
% 
% calctime = calc.t;
% ycalc = calc.y(1:ngene,:);
% 
% nexptime = length(exptime);
% ncalctime = length(calctime);
% 














% %==================
% %Step 1 of bilinear algorithm
% %==================
% 
% %Known
% V = diag(kmodel.U*kmodel.Unet);%Initial elementary flux
% X = kmodel.SpeciesConc;%Initial Concentrations
% S = kmodel.Mx;%Stoichiometric Matrix
% K = diag(kmodel.RateConst);

%Loop jreact <- 1:tnsteps
% react = log(X(S(:,jreact)<0));
% vj = ln(V(jreact));
% kj = vj - sum(react);
% K(jreact) = kj;
%end Loop

%Formulate Bi-Linear problem 

% kmodel.MConc;
% Vnet = kmodel.Unet;%Intial SS net reaction flux
% V = kmodel.U;%Initial elementary reaction flux (vj)
% X = kmodel.EnzSpecConc;%Initial Concentrations (C1, C2, etc)
% M = kmodel.M;%Stoichiometric matrix
% 
% elflux_v = zeros(tnsteps,1);
% for irxn = 1:tnsteps
%     metabs = X(M(:,irxn)<0);
%     
% end











% npoints = solverP.MaxDataPoints;
% nvars = size(Solution.initSS{2},1);
% 
% mean = zeros(nvars,1);
% mean(mean==0) = 1e-9;
% %Noise = randn(nvars,npoints).*1e-14;
% 
% NoisyData = Solution.initSS{2} + (randn(nvars,npoints).*1e-14);
   
 
%==========================================================================
%Plotting Dynamic data obtained from storeSolution.m
%==========================================================================
%Cycle through every pair of columns?
%Plot each case (initialSS, perturbation 1, etc) separately?

% genename = {'arcA'};
% pname = {};

% saveData = struct([]);
% figh = 0;
% %plotSSexpression(trnmodel,[],genename,pname,ngene,tnreg,saveData,allSolution)
% plotSSexpression(trnmodel,finalSSdata,genename,pname,ngene,tnreg,saveData);

% %GeneStruct.gname = {'arcA';'fnr'};
% GeneStruct = struct();
% GeneStruct.regulators = {{'ArcA';'FNR'};{'FNR';'OxygenRecp'}};
% GeneStruct.coeff = {[1e-7;1e-7];[1e-5;1e-5]};
% GeneStruct.srate = {[1e-3;1e-4];[1e-4;1e-4]};
% 
% [trnmodel] = changeparameter(trnmodel,defparval,regprotein,GeneStruct);

% genename = {'sdhA'};
% pname = {'OxygenRecp'};
% 
% %dbstop in savefile.m
% saveData.filename = 'Comparisons';
% plotSSexpression(trnmodel,[],genename,pname,ngene,tnreg,saveData,allSolution);

% linestyle = {'-';'-.';'--';':'};
% 
% %if exist('allSolution','var')
% all_sfnames = fieldnames(allSolution);
% nasf = length(all_sfnames);
% 
% n_gene = length(genename);
% for igene = 1:n_gene
%     tfg = strcmp(genename{igene},trnmodel.Gene);
%     if any(tfg)
%         
%         exp = sprintf('%s_(\\w+.?)',trnmodel.Gene{tfg});
%         celltemp = regexp(all_sfnames,exp,'tokens');
%         arraytemp = zeros(size(celltemp));
%         arraytemp(~cellfun('isempty',celltemp)) = 1;
%         gsfnames = all_sfnames(logical(arraytemp));
%         
%         kgsf = 1;
%         hcfig = zeros(length(gsfnames),1);
%         while kgsf <= length(gsfnames)  
%             hcfig(kgsf) =...
%             figure('Units','pixels','Position',[20 50 750 450],'Visible','on');
%             [nrows,ncols] = size(allSolution.(gsfnames{kgsf}));
%             
%             nexpts = ncols/2;
%             hline = zeros(ncols,1);
%             iexpt = 1;
%             while iexpt <= ncols               
%                 hline(iexpt) = ...
%                 line(allSolution.(gsfnames{kgsf})(:,iexpt),...
%                 allSolution.(gsfnames{kgsf})(:,iexpt+1)); 
%                 %Set Line Properties
%                 LProperty.LineStyle = '--';                
%                 LProperty.Color = rand(3,1)';
%                 setproperty(hcfig(kgsf),hline(iexpt),LProperty);
%                 
%                 
%                 
%             hold on
%             iexpt = iexpt + 2;
%             end 
%             hline(hline == 0) = [];
%             %Set Figure Properties
%             FProperty.Name = figname;
%             set(hcfig(kgsf),hcfig(kgsf),FProperty);
%             %Set Axes Properties
%             AProperty = struct();
%             set(hcfig(kgsf),gca(hcfig(kgsf)),AProperty);
%             
% 
%             kgsf = kgsf + 1;
%         end
%         
%     end
% end
 
    

%ODE solver call
%[time, dG] = ode23(@materialbalancefun,[tstart tend],initval);

   



