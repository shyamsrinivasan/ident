function plot2d_dynamic(xdata,ydata,axish,hf)
if nargin<4
    hf = figure;
end
if nargin<3
    axish = 0;
end

if axish
    axish.NextPlot = 'add';    
end


    