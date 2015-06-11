function [plotData,h_subfig] =...
printResults(trnmodel,ng,allSolution,conc,petconc,flux,printvar,varname)
%Plot Solution Curves
LineP = struct();
LineP.LineWidth = 2.0;
nexpt = length(fieldnames(allSolution));
% ColorSpec = setLineP(E{7});
for iexpt = 1:nexpt
    %Plot initial solution curve
    if iexpt == 1
        sim_name = 'init'; 
        LineP.Color = rgb('Green');
        LineP.DisplayName = sprintf('Initial');
        LineP.LineStyle = '--';
    elseif iexpt > 1
        %Select model
        sim_name = sprintf('sample_%d',iexpt-1); 
        %Plot solution curve
        LineP.Color = 'Black';
        LineP.LineStyle = '-';
        LineP.Displayname = sim_name;
    end
    if exist('h_subfig','var')
        [h_fig,h_subfig] =...
        plotSims(sim_name,trnmodel,ng,varname,allSolution,LineP,h_subfig);
    else
        [h_fig,h_subfig] =...
        plotSims(sim_name,trnmodel,ng,varname,allSolution,LineP);
    end        
end
setProperties(h_fig,h_subfig,allSolution.(sim_name));
%Steady state concentrations (all samples)
% var_n = trnmodel.Regulators;
% plotSSexpression(trnmodel,conc,var_n);
%Bifurcation Diagram
plotData = [];
% plotData = bifurcationPlot(ColorSpec,par_,nexpt,trnmodel,ng,varname,allfinalSS);
return
function [h_fig,h_subfig] =...
         plotSims(sim_name,trnmodel,ng,varname,allSolution,LineP,h_subfig)
%Plot solution curve                                                
    if isempty(findobj('type','figure','Name','Initial Solution'))
        h_fig = figure('Name','Initial Solution','Color',[1 1 1]);
    else
        h_fig = findobj('type','figure','Name','Initial Solution');
        figure(h_fig);
    end 
    if nargin > 6 || exist('h_subfig','var')
        [h_subfig] =...
        dynamicplot(trnmodel,ng,varname,allSolution.(sim_name),h_fig,LineP,h_subfig);
    else
        [h_subfig] =...
        dynamicplot(trnmodel,ng,varname,allSolution.(sim_name),h_fig,LineP);
    end
return
function setProperties(hfig,hsubfig,Solution,finalSS)
if nargin < 4
    %Solution call
    Xmin = 0;
    Xmax = max(Solution.t(:)/3600);
else
    %bifurcation call
    Xmin = min(finalSS.val(:,1));
    Xmax = max(finalSS.val(:,1));
end
if ~isempty(hsubfig)
    for ivar = 1:length(hsubfig)
%         Xmax = max(Solution.t(:)/3600);
%         Xdiv = sd_round(max(Solution.t(:))/10,1,1);
        %Common Axis Properties
        AxisP = struct();
        AxisP.XMinorTick = 'on';        
%         AxisP.XTick = 0:Xdiv:Xmax;
        AxisP.YMinorTick = 'on';
        AxisP.TickLength = [0.04,0.04];
        AxisP.XColor = [.1 .1 .1];
        AxisP.YColor = [.1 .1 .1];
        AxisP.XLim = [Xmin Xmax];
        AxisP.FontName = 'Lucida Sans';
        AxisP.FontSize = 12; 
        AxisP.FontWeight = 'bold';
        
        %Variable Specific Properties
        hca = findobj(hsubfig(ivar),'type','axes');
        set(get(hca,'XLabel'),'FontName','Lucida Sans',...
                              'FontSize',12,...
                              'FontWeight','bold');
%         set(get(hca,'XLabel'),);
        set(get(hca,'YLabel'),'FontName','Lucida Sans',...
                              'FontSize',12,...
                              'FontWeight','bold');   
%         set(get(hca,'YLabel'),);
        set(hfig,'CurrentAxes',hca);
        hline = findobj(hca,'type','line');
        Y = zeros(length(hline),2);
        iline = 1;
        while iline <= length(hline)
            Y(iline,1) = min(get(hline(iline),'YData'));
            Y(iline,2) = max(get(hline(iline),'YData'));
            iline = iline+1;
        end  
%         Ymin = sd_round(min(Y(:,1)),3,5);
%         Ymax = sd_round(max(Y(:,2)),3,5);
%         Ydiv = sd_round((max(Y(:,2))-min(Y(:,1)))/10,2,5);
%         AxisP.YLim = [Ymin Ymax];           
%         AxisP.YTick = Ymin:Ydiv:Ymax;
        if ~isempty(AxisP)
            set(hca,AxisP);
        end
    end
    legend(findobj(hsubfig(1),'type','axes'),'show');
    legend(findobj(hsubfig(1),'type','axes'),'boxoff');
end
