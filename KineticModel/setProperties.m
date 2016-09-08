function setproperties(plotype,hsfig,xlabel,ylabel,zlabel)
if nargin<5
    zlabel = {};
end
% set axis properties
set(get(gca,'YLabel'),'String',ylabel);  
set(get(gca,'YLabel'),'FontName','Arial');   
set(get(gca,'YLabel'),'FontSize',22); 
set(get(gca,'XLabel'),'String',xlabel);  
set(get(gca,'XLabel'),'FontName','Arial');   
set(get(gca,'XLabel'),'FontSize',22);
if plotype == 3
    set(get(gca,'ZLabel'),'String',zlabel);  
    set(get(gca,'ZLabel'),'FontName','Arial');   
    set(get(gca,'ZLabel'),'FontSize',22);
    axesP.ZColor = [.1 .1 .1];
end
axesP.FontName  = 'Arial';
axesP.FontSize = 22;
axesP.LineWidth = 1.5;
axesP.TickLength = [0.01 0.01];
axesP.XColor = [.1 .1 .1];
axesP.YColor = [.1 .1 .1];
if plotype == 3
    set(gca,axesP);
else
    set(hsfig,axesP);
end

% if ~isempty(hsubfig)
%     for ivar = 1:length(hsubfig)
%         Xmax = max(Solution.t(:));
%         Xdiv = sd_round(max(Solution.t(:))/10,1,1);
%         %Common Axis Properties
%         AxisP = struct();
% %         AxisP.XMinorTick = 'on';
%         AxisP.XTick = 0:Xdiv:Xmax;
% %         AxisP.YMinorTick = 'on';
%         AxisP.TickLength = [0.02,0.02];
%         AxisP.XColor = [.4 .4 .4];
%         AxisP.YColor = [.4 .4 .4];
%         AxisP.XLim = [0 Xmax];
%         AxisP.FontName = 'CMU Serif';
%         AxisP.FontSize = 24; 
% 
%         %Variable Specific Properties
%         if hsubfig(ivar)
%             hca = findobj(hsubfig(ivar),'type','axes');
%             set(hfig,'CurrentAxes',hca);
%             hline = findobj(hca,'type','line');
%             Y = zeros(length(hline),2);
%             iline = 1;
%             while iline <= length(hline)
%                 Y(iline,1) = min(get(hline(iline),'YData'));
%                 Y(iline,2) = max(get(hline(iline),'YData'));
%                 iline = iline+1;
%             end  
%     %         Ymin = min(Y(:,1));
%     %         Ymax = max(Y(:,2));
%     %         Ymin = sd_round(min(Y(:,1)),3,5);
%     %         Ymax = sd_round(max(Y(:,2)),3,5);
%     %         Ydiv = max(Y(:,2))-min(Y(:,1))/10;
%     %         Ydiv = sd_round((max(Y(:,2))-min(Y(:,1)))/10,2,5);
%     %         AxisP.YLim = [Ymin Ymax];           
%     %         AxisP.YTick = Ymin:Ydiv:Ymax;
%             if ~isempty(AxisP)
%                 set(hca,AxisP);
%             end
%         end
%     end
%     legend(findobj(hsubfig(1),'type','axes'),'show','Location','Best');
%     legend(findobj(hsubfig(1),'type','axes'),'boxoff');
% end 