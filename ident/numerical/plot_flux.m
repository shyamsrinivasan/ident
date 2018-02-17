function plot_flux(data, hfig, haxis, nplots, nrows, ncolumns)
for iplot = 1:nplots
    haxis(iplot) = subplot(nrows, ncolumns, iplot);
    xdata = data(iplot).t;
    ydata = data(iplot).flux;
    axis_labels = {'Time (s)','Flux a.u.'};
    plot2d(xdata,ydata,haxis(iplot),hfig,axis_labels);
end