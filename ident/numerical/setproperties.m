function setproperties(hsfig,axis_labels)

% set axis properties
if ~isempty(axis_labels)
    set(get(gca,'YLabel'),'String',axis_labels{2});  
    set(get(gca,'YLabel'),'FontName','Arial');   
    set(get(gca,'YLabel'),'FontSize',22); 
    set(get(gca,'XLabel'),'String',axis_labels{1});  
    set(get(gca,'XLabel'),'FontName','Arial');   
    set(get(gca,'XLabel'),'FontSize',22);
end

axesP.FontName  = 'Arial';
axesP.FontSize = 22;
axesP.LineWidth = 1.5;
axesP.TickLength = [0.01 0.01];
axesP.XColor = [.1 .1 .1];
axesP.YColor = [.1 .1 .1];
set(hsfig,axesP);
