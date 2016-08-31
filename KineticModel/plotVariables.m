function [hf,hsf] = plotVariables(xdata,ydata)
[nr,nc,np] = size(ydata);
if np>1
    for ip = 1:np
        subplot(2,2,1);
        plot(xdata,ydata(:,1,ip));
        ylabel('super Enzyme E');
        hold on
        subplot(2,2,2);
        plot(xdata,ydata(:,2,ip));
        ylabel('PEP');
        hold on
        subplot(2,2,3);
        plot(xdata,ydata(:,3,ip));
        ylabel('FBP');
        xlabel('time');
        hold on
    end
else
    subplot(2,2,1);
    plot(xdata,ydata(:,1));
    ylabel('super Enzyme E');    
    subplot(2,2,2);
    plot(xdata,ydata(:,2));
    ylabel('PEP');    
    subplot(2,2,3);
    plot(xdata,ydata(:,3));
    ylabel('FBP');
    xlabel('time');    
end
    