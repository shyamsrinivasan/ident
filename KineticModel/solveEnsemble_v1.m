%function [inSolution,enSolution] =...
%          solveEnsemble(ensb,model,variable,pertb,SolverOptions)
%solve an ensemble of models 
function [inSolution,enSolution] =...
         solveEnsemble(ensb,model,variable,pertb,SolverOptions)
nmodels = length(fieldnames(ensb));
ngraph(1) = 0;
saveData.dirname = 'C:\Users\shyam\Documents\Courses\CHE 1125 Project\Results\Kmodel';
%Initial Model Perturbation
fprintf('Initial Model Integration to obtain initial SS\n');
for isample = 1:nmodels
    mname = sprintf('model%d',isample);
    inSolution.(mname) = struct([]);
    fprintf('%s\n',mname);
%     saveData.filename = mname;
    %simulate models first to get initial SS
    Solution = callODEsolver(model,ensb.(mname),variable,inSolution.(mname),SolverOptions);
    inSolution.(mname) = Solution;   
%     savefile(Solution,mname,saveData);
    close all
end
%Enzyme Perturbations
fprintf('Model Perturbations\n');
npertb = length(fieldnames(pertb));
for isample = 1:nmodels
    mname = sprintf('model%d',isample);
    fprintf('%s\n',mname);
    for ipertb = 1:npertb
        pname = sprintf('pertb%d',ipertb);
        fprintf('Perturbation %d\n',ipertb);
        [Solution] = pertbEnzyme(pertb.(pname),model,ensb.(mname),...
                                 variable,inSolution.(mname),SolverOptions);
        enSolution.(pname).(mname) = Solution;
    end
    close all
end
%save all data
% [nvar,ntpts] = size(enSolution.pertb1.model1.y);
% save_data.t = enSolution.pertb1.model1.t;
% for isample = 1:nmodels
%     save_data.y = zeros(npertb*nvar,ntpts);
%     mname = sprintf('model%d',isample);
%     saveData.filename = mname;                    
%     ipertb = 0;
%     while ipertb < npertb
%         pname = sprintf('pertb%d',ipertb+1);
%         save_data.y(nvar*ipertb+1:nvar*(ipertb+1),:)=enSolution.(pname).(mname).y;
%         ipertb = ipertb+1;
%     end
%     save_data.y = [inSolution.(mname).y;save_data.y];
%     savefile(save_data,mname,saveData);
%     status = selectData_kmodel(save_data,model,mname);
% end

%Plots - Comparison between models
LineP = struct();
LineP.LineWidth = 1.5;
ColorSpec = setLineP(nmodels);
varname = {'M1','M2','M3'};
[hfig,hsubfig] = compareEnsemble(model,npertb,varname,enSolution,ColorSpec,LineP);
setProperties(hfig,hsubfig,enSolution.pertb1.model1);
%set default Line Properties

%Superceded by compareEnsemble
h_enfig = struct([]);
h_efig = zeros(npertb,1);
for isample = 1:nmodels
    %Select model
    mname = sprintf('model%d',isample);
%     h_enfig(1).(mname) = struct([]);
    %Plot initial solution curve
    LineP.Color = ColorSpec{isample};
    LineP.DisplayName = sprintf('model %d',isample);
    if isempty(findobj('type','figure','Name','Initial Solution'))
        h_ifig = figure('Name','Initial Solution');
    else
        figure(h_ifig);
    end      
    if exist('h_infig','var')
        [h_infig] = plotEnsemble(model,ngraph,varname,inSolution.(mname),h_ifig,LineP,h_infig);
    else
        [h_infig] = plotEnsemble(model,ngraph,varname,inSolution.(mname),h_ifig,LineP);
    end
    %Plot perturbation Solution
    for ipertb = 1:npertb
        pname = sprintf('pertb%d',ipertb);
        fname = sprintf('Perturbation Solution %s',pname);
        if isempty(findobj('type','figure','Name',fname))            
            h_efig(ipertb) = figure('Name',fname);
        else
            figure(h_efig(ipertb));
        end        
        
        if isfield(h_enfig,pname)  
            if isfield(h_enfig.(pname),mname)
                [h_enfig.(pname).(mname)] =...
                plotEnsemble(model,ngraph,varname,...
                             enSolution.(mname).(pname),...
                             h_efig(ipertb),LineP,h_enfig.(pname).(mname));
            else
                [h_enfig.(pname).(mname)] =...
                plotEnsemble(model,ngraph,varname,...
                             enSolution.(mname).(pname),...
                             h_efig(ipertb),LineP);
            end
        else
            [h_enfig(1).(pname).(mname)] =...
            plotEnsemble(model,ngraph,varname,...
                         enSolution.(mname).(pname),...
                         h_efig(ipertb),LineP);
        end
       
        
    end
end
%Not working below this point - changes needed - December 1st 2014
%Initial Solution Figure Properties
setProperties(h_ifig,h_infig,inSolution.(mname));
for ifig = 1:npertb
    pname = sprintf('pertb%d',ifig);
%     for isample = 1:nmodels
%         mname = sprintf('model%d',isample);
        setProperties(h_efig(ifig),h_enfig.(pname).(mname),enSolution.(mname).(pname));
%     end
end

%Plots - Comparison between initial & perturbation curves

        

for isample = 1:nmodels
    mname = sprintf('model%d',isample);
    figname = sprintf('Model %d Comparison',isample);        
    if isempty(findobj('type','figure','Name',figname))
        h_cfig = figure('Name',figname);
    else
        figure(h_cfig);
    end 
    LineP.Color = ColorSpec{1};
    LineP.DisplayName = sprintf('Initial');
    h_csubfig = plotEnsemble(model,ngraph,varname,inSolution.(mname),h_cfig,LineP);
    LineP.Color = ColorSpec{end};
    LineP.DisplayName = sprintf('Perturbation 1');
    h_csubfig = plotEnsemble(model,ngraph,varname,enSolution.(mname),h_cfig,LineP,h_csubfig);
    setProperties(h_cfig,h_csubfig,inSolution.(mname));
end
return

function ColorSpec = setLineP(nlines)
ColorSpec = cell(nlines,1);
for icolor = 1:nlines
    fac = icolor/(1+icolor);
    ColorSpec{icolor} = [1/2*fac 3/4*fac 1/3*fac];  
end

function setProperties(hfig,hsubfig,Solution)
if ~isempty(hsubfig)
    for ivar = 1:length(hsubfig)
        Xmax = max(Solution.t(:));
        Xdiv = sd_round(max(Solution.t(:))/10,1,1);
        %Common Axis Properties
        AxisP = struct();
%         AxisP.XMinorTick = 'on';
        AxisP.XTick = 0:Xdiv:Xmax;
%         AxisP.YMinorTick = 'on';
        AxisP.TickLength = [0.02,0.02];
        AxisP.XColor = [.2 .2 .2];
        AxisP.YColor = [.2 .2 .2];
        AxisP.XLim = [0 Xmax];
        AxisP.FontName = 'Courier';
        AxisP.FontSize = 12; 

        %Variable Specific Properties
        hca = findobj(hsubfig(ivar),'type','axes');
        set(hfig,'CurrentAxes',hca);
        hline = findobj(hca,'type','line');
        Y = zeros(length(hline),2);
        iline = 1;
        while iline <= length(hline)
            Y(iline,1) = min(get(hline(iline),'YData'));
            Y(iline,2) = max(get(hline(iline),'YData'));
            iline = iline+1;
        end  
        Ymin = sd_round(min(Y(:,1)),3,5);
        Ymax = sd_round(max(Y(:,2)),3,5);
%         Ydiv = max(Y(:,2))-min(Y(:,1))/10;
%         Ydiv = sd_round((max(Y(:,2))-min(Y(:,1)))/10,2,5);
        AxisP.YLim = [Ymin Ymax];           
%         AxisP.YTick = Ymin:Ydiv:Ymax;
        if ~isempty(AxisP)
            set(hca,AxisP);
        end
    end
    legend(findobj(hsubfig(1),'type','axes'),'show','Location','Best');
    legend(findobj(hsubfig(1),'type','axes'),'boxoff');
end

%choose selected solution subplot
%plot solution



    
    
    %plot each solution curve (for one or more metabolites)?
    
%Plot enzyme perturbation results
%Save file information
% saveInfo.dirname = 'C:\Users\shyam\Documents\Courses\CHE 1125 Project\Results\Kinetic_Model';
% saveInfo.filename = '';%sprintf('ExptCondition_%d',exptnum);

%figfname = sprintf('%s_',saveData.filename);
%savefile(Sol,[],saveInfo,hfig,figfname);
% close all