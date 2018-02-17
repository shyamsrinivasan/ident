function hline = plot2d(xdata, ydata, haxis, hfig, axis_labels)
if nargin<5
    axis_labels = [];
end
if nargin<4
    hfig = figure;    
end
set(0,'CurrentFigure',hfig);
set(hfig,'CurrentAxes',haxis);
set(haxis,'NextPlot','add');
hline = line(xdata,ydata);

setproperties(haxis,axis_labels);




