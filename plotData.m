function [hobj] = plotData(t,y,h,LineP)
if nargin<4
    LineP.Color = [];
    LineP.DisplayName = {};
    LineP.LineWidth = [];
    LineP.YLabel = {};
else
    if ~isfield(LineP,'Color')
        LineP.Color = [];
    end
    if ~isfield(LineP,'DisplayName')
        LineP.DisplayName = {};
    end
    if ~isfield(LineP,'LineWidth')
        LineP.LineWidth = 0.5;
    end
    if ~isfield(LineP,'YLabel')
        LineP.YLabel = {};
    end
end
if length(h)>1
    hfig = h(1);
    hsf = h(2);
else
    hfig = h(1);
    hsf = [];
end

%axis labels
gca;
xlabel('Time');
% set(hsf,'XLabel','time');
if ~isempty(LineP.YLabel);
%     set(gca,'YLabel',LineP.YLabel);
    ylabel(LineP.YLabel);
end

%find axis corresponding to hsf
if ~isempty(hsf)
    hca = findobj(hsf,'type','axes');
else
    hca = gca;
end
set(hfig,'CurrentAxes',hca);

%plot data as a line to get object handle
hobj = line('XData',t,'YData',y);

%line color
if ~isempty(LineP.Color)
    set(hobj,LineP.Color);
else
    set(hobj,'Color',[0 .5 0]);
end

%line legend
if ~isempty(LineP.DisplayName)
    set(hobj,'DisplayName',LineP.DisplayName);
end

%line width
if ~isempty(LineP.LineWidth)
    set(hobj,'LineWidth',LineP.LineWidth);
else
    set(hobj,'LineWidth',0.5);
end


%Make sure more plots can be added
set(hca,'NextPlot','add');     
set(hsf,'NextPlot','add');

    



    

